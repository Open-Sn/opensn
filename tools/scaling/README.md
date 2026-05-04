# Node-to-Node Scaling Study

This is a node-to-node scaling study for OpenSn. For now, tets-only unstructured meshes with AAH sweep
are supported.

## 1. Generate inputs

``generate_scaling_study.py`` generates the meshes and launch scripts necessary for strong an weak scaling studies.

Usage:
```
python3 generate_scaling_study.py
  --type {strong,weak}  Type of scaling test to generate files for.
  --processor {cpu,gpu} Processor target for the generated study files.
                        Defaults to CPU.
  --engine {slurm,lc-flux,alcf-pbs}
                        Job submitting system. Defaults to slurm.
```

Examples:
```
python3 generate_scaling_study.py --type=strong --processor=cpu --engine=slurm
```
```
python3 generate_scaling_study.py --type=weak --processor=gpu --engine=lc-flux
```

After running the script, a folder ``output/{strong/weak}_{cpu/gpu}`` will appear.
Change directory into that folder for the next step.

## 2. Submitting jobs

Inside of the output folder (``output/*/``), execute:
```
source submit_jobs.sh
```
and wait for all the jobs to finish.

When all jobs have been finished, go back the ``scaling`` folder (i.e. exit the ``output/`` folder).

## 3. Plot

To plot strong/weak scaling, run:
```
python3 plot_strong.py --dir=output/strong_cpu
```
```
python3 plot_weak.py --dir=output/weak_gpu
```

Result of the current run can be compared with or recorded with:
```
python plot_strong.py
  --history {none,comp,save}
                        History mode for the plot:
                            none (default, only plot current data),
                            comp (compare with history in the same plot without saving),
                            save (plot and overwrite current history value).
```
