MkDocs 
==================================

MkDocs is a fast, simple and downright gorgeous static site generator that's geared towards building project documentation. Documentation source files are written in Markdown, and configured with a single YAML configuration file. Start by reading the introductory tutorial, then check the User Guide for more information.

#. `MkDocs Official Website <https://www.mkdocs.org/>`_
#. `How to use CURL on Windows | How to test API with CURL | CURL Basics Step by Step <https://www.youtube.com/watch?v=8f9DfgRGOBo/>`_
#. `[Tutorial] - How to Create A GitHub Repo using the GitHub REST API | 2022 <https://www.youtube.com/watch?v=5djgwx9aWrg/>`_


install mkdocs
::

  pip install mkdocs
  
mkdocs new mkdocs-site
::

  mkdocs new mkdocs-site
  results:
  PS D:\github> mkdocs new mkdocs-site
  INFO    -  Creating project directory: mkdocs-site
  INFO    -  Writing config file: mkdocs-site\mkdocs.yml
  INFO    -  Writing initial docs: mkdocs-site\docs\index.md 

mkdocs serve
::

  mkdocs serve
  
pip install mkdocs-material
::

  pip install mkdocs-material

::

  $ git init
  $ git add .
  $ git commit -m "init"  
  
::
  
  git remote add origin git@github.com:eric2003/mkdocs-site.git
  git branch -M main
  git push -u origin main
  
::

  $ mkdir .github
  $ cd .github
  $ mkdir workflows
  $ cd workflows
  $ vim PublishMySite.yml  
  
PublishMySite.yml
::

  name: publish site
  on: # 在什么时候触发工作流
    push: # 在从本地main分支被push到GitHub仓库时
      branches:
        - main
    pull_request: # 在main分支合并别人提的pr时
      branches:
        - main
  jobs: # 工作流的具体内容
    deploy:
      runs-on: ubuntu-latest # 创建一个新的云端虚拟机 使用最新Ubuntu系统
      steps:
        - uses: actions/checkout@v2 # 先checkout到main分支
        - uses: actions/setup-python@v2 # 再安装Python3和相关环境
          with:
            python-version: 3.x
        - run: pip install mkdocs-material # 使用pip包管理工具安装mkdocs-material
        - run: mkdocs gh-deploy --force # 使用mkdocs-material部署gh-pages分支  
  
github api access token curl
::

  Log in to your GitHub account.
  Click on your profile icon at the top right corner and select "Settings".
  Choose the "Developer settings" tab.
  Click on "Personal access tokens".
  Click on the "Generate new token" button.
  Enter your password to verify your identity.
  In the "Note" field, enter a brief description of the purpose of your access token.
  Select the permissions you want this token to have. For example, if you want to access your public repositories, you need to select the "public_repo" permission.
  Click the "Generate token" button." 

  
  
