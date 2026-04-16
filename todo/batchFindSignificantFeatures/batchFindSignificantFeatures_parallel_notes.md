# `batchFindSignificantFeatures` 并行行为与修改建议

文件：`DROMA_R/R/FuncPairBatchFeature.R`  
目标函数：`batchFindSignificantFeatures()`

## 结论

当前实现可以并行，但不能保证运行全程始终维持 `cores` 个活跃 worker。

原因不是 `future::plan(..., workers = cores)` 没生效，而是实际并行粒度是 `chunks`，不是单个 feature。因此真实并行上限是：

`min(cores, length(chunks))`

所以当：

- `feature2_list` 较短
- `chunk_size` 偏大
- 某些 chunk 提前结束

时，后台看到的活跃进程数少于设定的 `cores` 是正常现象。

## 当前实现定位

### 1. `cores` 确实传给了 future backend

函数中这段代码会根据平台设置并行 backend：

- Unix/Mac: `future::multicore`
- Windows: `future::multisession`

核心代码：

```r
future::plan(future::multicore, workers = cores)
```

或

```r
future::plan(future::multisession, workers = cores)
```

这说明 `cores` 参数本身确实参与了并行配置。

### 2. 真正并行的是 chunk，不是 feature

后续逻辑先按 `chunk_size` 把 `feature2_list` 切块：

```r
chunks <- split(seq_along(feature2_list),
               ceiling(seq_along(feature2_list) / chunk_size))
```

然后对 `chunks` 做：

```r
furrr::future_map(chunks, function(indices) {
  lapply(indices, worker_function)
})
```

这意味着每个 future 处理的是一个 chunk，chunk 内部仍然是串行 `lapply()`。

因此：

- 如果 `length(chunks) < cores`，就不可能跑满 `cores`
- 即使 `length(chunks) >= cores`，也可能因为 chunk 完成时间不同，导致后半程活跃 worker 下降

## 为什么会出现“设置 8 核，但后台只看到 2-4 个”

### 情况 1：feature 数不够多

当 `n_features < 100` 时：

```r
chunk_size <- max(10, ceiling(n_features / cores))
```

例如：

- `n_features = 30`
- `cores = 8`

则：

- `chunk_size = 10`
- `length(chunks) = 3`

所以最多只会有 3 个并行 worker。

### 情况 2：chunk 太粗

当每个 chunk 很大时，任务分配粒度粗，worker 利用率会下降：

- 有些 worker 很早做完
- 剩下少数 worker 处理“重 chunk”
- 监控上就会看到活跃 worker 数变少

### 情况 3：并行的是进程，不一定是你看到的“线程”

在 Unix/Mac 上，当前用的是 `future::multicore`，这是 fork 模式，更接近多进程而不是单进程内多线程。

所以系统监控里如果你看的口径是 thread 数，不一定和 `cores` 一一对应。

### 情况 4：未 preload 时可能有 I/O/数据库争用

代码里已经有提示：

```r
Note: Parallel processing without preloaded data may have database contention.
```

如果 `preloaded_feas2` 为空，各 worker 可能争抢底层数据读取资源。这样即便 worker 数不少，CPU 也不一定打满。

## 修改建议

### 建议 1：增加并行诊断日志

建议在切块完成后打印以下信息：

- `n_features`
- `cores`
- `chunk_size`
- `n_chunks`
- `actual_parallelism = min(cores, n_chunks)`

建议日志示例：

```r
n_chunks <- length(chunks)
actual_parallelism <- min(cores, n_chunks)
message(sprintf(
  "Parallel setup: n_features=%d, cores=%d, chunk_size=%d, n_chunks=%d, max_active_workers=%d",
  n_features, cores, chunk_size, n_chunks, actual_parallelism
))
```

价值：

- 能直接解释为什么没跑满 `cores`
- 方便判断问题出在 chunking 还是 backend

### 建议 2：减小小任务场景下的 chunk_size

当前小任务规则：

```r
if (n_features < 100) {
  max(10, ceiling(n_features / cores))
}
```

这个 `max(10, ...)` 会导致小任务时 chunk 过粗。

可考虑改成更细粒度，例如：

```r
if (n_features < 100) {
  max(1, ceiling(n_features / (cores * 2)))
}
```

或者更保守一点：

```r
if (n_features < 100) {
  max(2, ceiling(n_features / (cores * 2)))
}
```

效果：

- 小规模任务更容易接近 `cores`
- 调度更均匀

代价：

- future 调度开销会上升
- 但对“明明设了多核却看起来没并行”的问题更友好

### 建议 3：限制 chunk_size 不要让 `n_chunks < cores`

如果目标是“尽量打满 worker”，可以在生成 chunk 前增加约束：

```r
chunk_size <- min(chunk_size, ceiling(n_features / cores))
chunk_size <- max(1, chunk_size)
```

这样至少能保证：

```r
length(chunks) >= cores
```

前提是 `n_features >= cores`。

注意：

- 这会偏向“打满 worker”
- 不一定是总耗时最优
- 但更符合用户对 `cores` 的直觉预期

### 建议 4：提供显式参数控制 chunk 策略

建议新增可选参数，例如：

- `chunk_size = NULL`
- 或 `tasks_per_core = 4`

示例：

```r
batchFindSignificantFeatures(..., cores = 8, tasks_per_core = 8)
```

然后内部计算：

```r
chunk_size <- max(1, ceiling(n_features / (cores * tasks_per_core)))
```

价值：

- 用户可根据数据规模调节并行粒度
- 比硬编码三段式规则更灵活

### 建议 5：如果重点是吞吐量，优先配合 `preloaded = TRUE`

如果 `feature2_type` 支持 preload，建议在文档和 message 中更明确提示：

- `cores > 1` 时优先使用 `preloaded = TRUE`
- 否则并行收益可能被 I/O/数据库访问抵消

这不是“并行数量变少”的直接原因，但会影响你观察到的实际 CPU 利用率。

## 推荐改法

如果目标是：

“让并行行为更符合 `cores` 的直觉，同时尽量少改现有结构”

建议采用下面这组最小修改：

1. 保留 `future::multicore` / `multisession` 逻辑不变
2. 保留 `furrr::future_map(chunks, ...)` 结构不变
3. 增加并行诊断日志
4. 调整 `n_features < 100` 时的 chunk 规则，让 chunk 更细
5. 增加一个保护，避免 `n_chunks` 明显小于 `cores`

示意实现：

```r
if (n_features < 100) {
  chunk_size <- max(1, ceiling(n_features / (cores * 2)))
} else if (n_features < 1000) {
  chunk_size <- ceiling(n_features / (cores * 4))
} else {
  chunk_size <- max(50, min(200, ceiling(n_features / (cores * 8))))
}

chunk_size <- min(chunk_size, ceiling(n_features / cores))
chunk_size <- max(1, chunk_size)

chunks <- split(seq_along(feature2_list),
               ceiling(seq_along(feature2_list) / chunk_size))

n_chunks <- length(chunks)
actual_parallelism <- min(cores, n_chunks)
message(sprintf(
  "Parallel setup: n_features=%d, cores=%d, chunk_size=%d, n_chunks=%d, max_active_workers=%d",
  n_features, cores, chunk_size, n_chunks, actual_parallelism
))
```

## 验证建议

修改后建议至少验证以下场景：

### 场景 1：小任务

- `n_features = 20~50`
- `cores = 8`

观察：

- 日志中的 `n_chunks`
- 系统监控中的活跃 worker 数

### 场景 2：中等任务

- `n_features = 200~500`
- `cores = 8`

观察：

- 是否比旧版本更稳定地维持多 worker 活跃
- 总运行时间是否明显变差

### 场景 3：大任务

- `n_features > 2000`
- 比较 `preloaded = TRUE/FALSE`

观察：

- CPU 利用率
- 总耗时
- I/O 等待是否明显

## 一句话总结

当前代码“能并行”，但并行单位是 chunk，不是 feature，所以实际活跃 worker 数常常会低于 `cores`。如果想让运行表现更接近“我设了几核就尽量跑几核”，应优先调整 chunk 策略，并增加并行诊断日志。
