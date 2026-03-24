# Node-to-Node Scaling Study

This is a node-to-node scaling study for OpenSn. For now, tets-only unstructured meshes are
supported.

GPU-accelerated scaling study can be performed by setting `ngpus` in `generate_scaling_study.py`.
This input will overwrite the number of tasks per nodes and cores per task.

## Strong scaling

To generate inputs for strong scaling study:

```
python3 generate_scaling_study --type=strong
```

Jobs can be submitted with:

```
cd output/strong
source submit_jobs.sh
```

Once results are ready:
```
python3 plot_strong.py output/strong
```

## Weak scaling

To generate inputs for weak scaling study:

```
python3 generate_scaling_study --type=weak
```

Jobs can be submitted with:

```
cd output/weak
source submit_jobs.sh
```

Once results are ready:
```
python3 plot_weak.py output/weak
```
