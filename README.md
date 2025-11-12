in pycharm select new project, and select conda NOT pip.
```commandline
conda install dftd3-python -c conda-forge
```

create a folder (say "water"), write code. then use terminal to run 
```commandline
cd water
python dft_water.py
```

this way, the timer.dat, water.molden, and water_optimized.molden will be inside the same directory, 
hence they won't get mixed.
