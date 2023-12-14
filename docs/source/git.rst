Git
==================================

Git
---------------------------------
#. `Git 教程 <https://www.runoob.com/git/git-push.html>`_
#. `Git系列链接整理 <https://zhuanlan.zhihu.com/p/503890935/>`_
#. `git常用命令 <http://www.pedestrian.com.cn/misc/git/git_normal.html#id2>`_
 

::

  git --version
  Typical results:
  PS C:\Users\eric> git --version
  git version 2.43.0.windows.1
  
查看远程仓库
::

  git remote -v
  Typical results:
  PS D:\github\OneFLOW> git remote -v
  origin  https://github.com/eric2003/OneFLOW.git (fetch)
  origin  https://github.com/eric2003/OneFLOW.git (push)  

::

  git add .
  git commit -m "first update"
  git push origin main
  
::

  git config --list --show-origin  
  Typical results:
  file:C:/Users/eric/.gitconfig
  
#. `git clone with SSH keys 或 HTTPS 設定步驟  <https://tsengbatty.medium.com/git-%E8%B8%A9%E5%9D%91%E7%B4%80%E9%8C%84-%E4%BA%8C-git-clone-with-ssh-keys-%E6%88%96-https-%E8%A8%AD%E5%AE%9A%E6%AD%A5%E9%A9%9F-bdb721bd7cf2/>`_

Check for an existing SSH key
::

  ls -al ~/.ssh
  Typical results:
  $ ls -al ~/.ssh
  ls: cannot access '/c/Users/eric/.ssh': No such file or directory
  
If you don't see any output or that directory doesn't exist (you get a No such file or directory message), then run:
::

  mkdir $HOME/.ssh  
  
  
Then generate a new set of keys with:
::

  ssh-keygen -t rsa -b 4096 -C your@email.com
  
#. `如何清空ssh私钥对应的passphrase？ <https://mingda.dev/2022/05/23/removing-passphrase-from-ssh-private-key/>`_
::

  ssh-keygen -p

git clone REPOSITORY.git by ssh
::

  git clone git@github.com:USERNAME/REPOSITORY.git

  
The authenticity of host 'github.com (192.30.255.113)' can't be established.
::

  ssh-keyscan github.com >> ~/.ssh/known_hosts
  
測試是否設定成功 & 設定 SSH agent
::

  ssh -T git@github.com
  PS C:\Users\eric> ssh -T git@github.com
  The authenticity of host 'github.com (192.30.255.113)' can't be established.
  ED25519 key fingerprint is ****.
  This key is not known by any other names
  Are you sure you want to continue connecting (yes/no/[fingerprint])? yes
  Warning: Permanently added 'github.com' (ED25519) to the list of known hosts.
  Hi eric2003! You've successfully authenticated, but GitHub does not provide shell access.
  
  
Recv failure: Connection was reset 
::

  git config --global --unset http.proxy
  git config --global --unset https.proxy
  ipconfig /flushdns

::
  
  git config --global http.proxy http://127.0.0.1:1234
  git config --global https.proxy http://127.0.0.1:1234


Download your recovery codes
You can use recovery codes as a second factor to authenticate in case you lose access to your device. We recommend saving them with a secure password manager such as 1Password, Authy, or Keeper.


git tag
::

  git tag 0.0.4
  
git branch
::

  git branch

git log
::

  git log
  Typical results:
  PS D:\work\github_work\Foo> git log
  commit 9b79604a930bf3f02f4eb9647b4316257d31866c (HEAD -> main, origin/main, origin/HEAD)
  Author: eric <****@hotmail.com>
  Date:   Tue Dec 5 18:55:47 2023 +0800
  
      add Foo:Foo ALIAS
  
  commit 0c34e18006d72918b170d91b4508690ae94f7517 (tag: 0.0.3)
  Author: eric <****@hotmail.com>
  Date:   Mon Dec 4 22:17:12 2023 +0800
  
      add config.cmake.in
  
  commit 7fa1ed87c1c8efc8f95d67cbcf7574731270162b (tag: 0.0.2)
  Author: eric <****@hotmail.com>
  Date:   Tue Nov 28 18:27:21 2023 +0800
  
      add install
  
  commit c3f207cc80128b37acfa93a51731a0515a0bb5ad (tag: 0.0.1)
  Author: eric <****@hotmail.com>
  Date:   Tue Nov 28 15:59:27 2023 +0800
  
      add foo
  
  commit 8e5d2f82edd4c875076094360ca8d710b909f38a
  Author: eric <****@hotmail.com>
  Date:   Tue Nov 28 15:28:46 2023 +0800
  
      Initial commit  

quit git log
::

  q
  
git log pretty
::

  git log --pretty=oneline --abbrev-commit
  Typical results:
  PS D:\work\github_work\Foo> git log --pretty=oneline --abbrev-commit
  9b79604 (HEAD -> main, origin/main, origin/HEAD) add Foo:Foo ALIAS
  0c34e18 (tag: 0.0.3) add config.cmake.in
  7fa1ed8 (tag: 0.0.2) add install
  c3f207c (tag: 0.0.1) add foo
  8e5d2f8 Initial commit
  
git tag
::

  git tag -a v0.0.4 9b79604 
  Typical results:
  PS D:\work\github_work\Foo> git tag -a v0.0.4 9b79604
  PS D:\work\github_work\Foo> git tag
  0.0.1
  0.0.2
  0.0.3
  v0.0.4
  
git tag quit
::

  保存并退出：
 （1）按 **Esc**键退出编辑模式，英文模式下输入 :wq ，然后回车(write and quit)。
 （2）按 Esc 键退出编辑模式，大写英文模式下输入 ZZ ，然后回车。

  不保存退出：
  按 **Esc**键退出编辑模式，英文模式下输入 :q! ，然后回车。
  按 **Esc**键退出编辑模式，英文模式下输入 :qa! ，然后回车。
  
git tag modify
::

  git tag new-tag-name old-tag-name
  git tag -d old-tag-name 

git tag modify example
::

  git tag 0.0.4 v0.0.4
  git tag -d v0.0.4
  Typical results:
  PS D:\work\github_work\Foo> git tag 0.0.4 v0.0.4
  PS D:\work\github_work\Foo> git tag -d v0.0.4
  Deleted tag 'v0.0.4' (was 623a70e)
  PS D:\work\github_work\Foo> git tag
  0.0.1
  0.0.2
  0.0.3
  0.0.4
  
git remote -v
::

  PS D:\work\github_work\Foo> git remote -v
  origin  git@github.com:eric2003/Foo.git (fetch)
  origin  git@github.com:eric2003/Foo.git (push)  
  
git push --tags
::

  git push --tags
  Typical results:
  PS D:\work\github_work\Foo> git push --tags
  Enumerating objects: 1, done.
  Counting objects: 100% (1/1), done.
  Writing objects: 100% (1/1), 162 bytes | 162.00 KiB/s, done.
  Total 1 (delta 0), reused 0 (delta 0), pack-reused 0
  To github.com:eric2003/Foo.git
   * [new tag]         0.0.4 -> 0.0.4  
   
gist.github.com
::

  https://gist.github.com/eric2003
  
example
::

   git add .
   git commit -m "add version"
   git tag -a v0.0.5 -m "New release for v0.0.5"
   git push origin main
   git push --tags