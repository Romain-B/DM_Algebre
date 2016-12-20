# DM_Algebre

##Launching an IPython notebook 

Check that you have [Jupyter](http://jupyter.org/) installed by typing `jupyter --version` in a terminal

If not, you can install it with `pip` :
```
$ sudo pip install --upgrade pip
$ suso pip install jupyter
```
Then, open a new terminal (it's going to be locked afterwards until you close the notebook server)
and type 
```
$ sudo jupyter notebook
```
*note : the `sudo` is a prevention measure since sometimes, you need root authorization to create files (such as ipython notebooks).*

Once you've done so, a server log should start rolling out to which you don't really need top pay attention to, other than the exit instructions which are basically `Ctrl+C` and `y` to confirm.

Your browser should open up with a navigation window to browse through your files, just go to the notebook and launch it.
If your browser doesn't launch the jupyter interface, got to this adress : [http://localhost:8888/](http://localhost:8888/). 
