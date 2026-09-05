Python
==================================

Online Python Code Editor
---------------------------------
#. `Online Python Code Editor <https://pynative.com/online-python-code-editor-to-execute-python-code/>`_

Python Tutorial
---------------------------------
#. `Python Tutorial <https://www.w3schools.com/python/>`_
#. `Python for Scientific Computing <https://sbu-python-class.github.io/python-science/Introduction.html>`_


uninstall pip:
::

   python -m pip uninstall pip
   py.exe -m pip uninstall pip
   py.exe -m ensurepip –upgrade

install pip:
::

  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
  py.exe get-pip.py
  


ModuleNotFoundError: No module named 'scipy':
::

    pip install scipy
	
ModuleNotFoundError: No module named 'pyamg':
::

    pip install pyamg	
	
ModuleNotFoundError: No module named 'tqdm':
::

    pip install tqdm

How to check python installation path   
::

    import sys
    print(sys.exec_prefix)

::

    py -0p

Typical output:
::

    PS D:\work\python_work> py.exe -0p
      -V:3.11 *        C:\Users\eric\AppData\Local\Programs\Python\Python311\python.exe
      -V:ContinuumAnalytics/Anaconda39-64 C:\Users\eric\miniconda3\python.exe

Check Installed Modules in Python

::

    py.exe -m pip list
	
Install the classic Jupyter Notebook with:

::

    pip install jupyter

.. toctree::
   :maxdepth: 1

   /python/matplotlib
   
::

  import os
  import shutil
  
  # 获取当前目录
  current_directory = os.getcwd()
  
  # 遍历当前目录下的所有子目录
  for root, dirs, files in os.walk(current_directory):
      for directory in dirs:
          if directory == "build":
              build_directory = os.path.join(root, directory)
              print(f"Deleting contents of {build_directory}")
              # 删除子目录中的所有内容
              shutil.rmtree(build_directory)
  
  print("Deletion complete")

::
  
  import os
  import shutil
  
  def remove_git_folder(path):
      os.system(f'rd /s /q {path}')
  
  def delete_build_directories(root_dir):
      print(f'root_dir {root_dir}')
      for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
          for dirname in dirnames:
              if dirname == 'build':
                  print(f'dirpath {dirpath} dirname {dirname}')
                  build_dir = os.path.join(dirpath, dirname)
                  print(f'build_dir {build_dir}')
                  try:
                      shutil.rmtree(build_dir)
                      print(f'Deleted {build_dir}')
                  except OSError as e:
                      if e.errno == 13:
                          print(f'Permission denied: {build_dir}')
                          remove_git_folder(build_dir)
                      elif e.errno == 16:
                          print(f'Directory not empty: {build_dir}')
                      else:
                          print(f'Error deleting {build_dir}: {e}')                
  
  # 获取当前目录
  current_directory = os.getcwd()
  modernCMakeDir = r'd:\work\modern_cmake_work\ModernCMake\codes'
  
  print(f"current_directory= {current_directory}")
  print(f"modernCMakeDir= {modernCMakeDir}")
  
  # 调用函数并传入要删除目录的根目录
  delete_build_directories(modernCMakeDir)
  input()
