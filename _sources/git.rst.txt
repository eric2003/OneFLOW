Git
==================================

Git
---------------------------------
#. `Git 教程 <https://www.runoob.com/git/git-push.html>`_
#. `Git系列链接整理 <https://zhuanlan.zhihu.com/p/503890935/>`_
#. `git常用命令 <http://www.pedestrian.com.cn/misc/git/git_normal.html#id2>`_
#. `github actions 从入门到精通系列链接整理 <https://zhuanlan.zhihu.com/p/388642124/>`_
 

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
   
new-branch
::

  git checkout -b new-branch
  git checkout -b develop
   
GitHub Actions
--------------------------------------------
#. `Using GitHub Actions with C++ and CMake <https://cristianadam.eu/20191222/using-github-actions-with-c-plus-plus-and-cmake/>`_
#. `GitHub Actions Step by Step DEMO for Beginners <https://www.youtube.com/watch?v=ylEy4eLdhFs/>`_
#. `GitHub Actions Tutorial | From Zero to Hero in 90 minutes (Environments, Secrets, Runners, etc) <https://www.youtube.com/watch?v=TLB5MY9BBa4/>`_
#. `GitHub Actions: Approvals, Environments and Visualization DEEP DIVE <https://www.youtube.com/watch?v=w_37LDOy4sI/>`_
#. `GitHub Actions Secrets: Security Best Practices <https://www.youtube.com/watch?v=2yHRq7aWDKM/>`_
#. `使用 GitHub Actions 实现持续集成与持续部署 <https://www.bilibili.com/video/BV1pf4y1776s/>`_


Contexts
============
#. `Contexts <https://docs.github.com/en/actions/learn-github-actions/contexts/>`_

echo
::

  ubuntu:
  echo "GITHUB_ACTION=$GITHUB_ACTION."
  windows:
  echo "GITHUB_ACTION=$env:GITHUB_ACTION."


Github Project Info
::

  name: GitHub Default environment variables
  on: [push]
  jobs:
    Display-GitHub-Default-environment-variables:
      runs-on: ubuntu-latest
      steps:
        - name: Default-environment-variables
          run: |    
            echo "CI=$CI."
            echo "RUNNER_OS=$RUNNER_OS."
            echo "GITHUB_ACTION=$GITHUB_ACTION."
            echo "GITHUB_ACTION_PATH=$GITHUB_ACTION_PATH."
            echo "GITHUB_ACTION_REPOSITORY=$GITHUB_ACTION_REPOSITORY."
            echo "GITHUB_ACTIONS=$GITHUB_ACTIONS."
            echo "GITHUB_ACTOR=$GITHUB_ACTOR."
            echo "GITHUB_ACTOR_ID=$GITHUB_ACTOR_ID."
            echo "GITHUB_API_URL=$GITHUB_API_URL."
            echo "GITHUB_BASE_REF=$GITHUB_BASE_REF."
            echo "GITHUB_ENV=$GITHUB_ENV."
            echo "GITHUB_EVENT_NAME=$GITHUB_EVENT_NAME."
            echo "GITHUB_EVENT_PATH=$GITHUB_EVENT_PATH."
            echo "GITHUB_GRAPHQL_URL=$GITHUB_GRAPHQL_URL."
            echo "GITHUB_HEAD_REF=$GITHUB_HEAD_REF."
            echo "GITHUB_JOB=$GITHUB_JOB."
            echo "GITHUB_OUTPUT=$GITHUB_OUTPUT."
            echo "GITHUB_PATH=$GITHUB_PATH."
            echo "GITHUB_REF=$GITHUB_REF."
            echo "GITHUB_REF_NAME=$GITHUB_REF_NAME."
            echo "GITHUB_REF_PROTECTED=$GITHUB_REF_PROTECTED."
            echo "GITHUB_REF_TYPE=$GITHUB_REF_TYPE."
            echo "GITHUB_REPOSITORY=$GITHUB_REPOSITORY."
            echo "GITHUB_REPOSITORY_ID=$GITHUB_REPOSITORY_ID."
            echo "GITHUB_REPOSITORY_OWNER=$GITHUB_REPOSITORY_OWNER."
            echo "GITHUB_REPOSITORY_OWNER_ID=$GITHUB_REPOSITORY_OWNER_ID."
            echo "GITHUB_RETENTION_DAYS=$GITHUB_RETENTION_DAYS."
            echo "GITHUB_RUN_ATTEMPT=$GITHUB_RUN_ATTEMPT."
            echo "GITHUB_RUN_ID=$GITHUB_RUN_ID."
            echo "GITHUB_RUN_NUMBER=$GITHUB_RUN_NUMBER."
            echo "GITHUB_SERVER_URL=$GITHUB_SERVER_URL."
            echo "GITHUB_SHA=$GITHUB_SHA."
            echo "GITHUB_STEP_SUMMARY=$GITHUB_STEP_SUMMARY."
            echo "GITHUB_TRIGGERING_ACTOR=$GITHUB_TRIGGERING_ACTOR."
            echo "GITHUB_WORKFLOW=$GITHUB_WORKFLOW."
            echo "GITHUB_WORKFLOW_REF=$GITHUB_WORKFLOW_REF."
            echo "GITHUB_WORKFLOW_SHA=$GITHUB_WORKFLOW_SHA."
            echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE."
            echo "RUNNER_ARCH=$RUNNER_ARCH."
            echo "RUNNER_DEBUG=$RUNNER_DEBUG."
            echo "RUNNER_NAME=$RUNNER_NAME."
            echo "RUNNER_OS=$RUNNER_OS."
            echo "RUNNER_TEMP=$RUNNER_TEMP."
            echo "RUNNER_TOOL_CACHE=$RUNNER_TOOL_CACHE."

Results
::

  Run echo "CI=$CI."
  CI=true.
  RUNNER_OS=Linux.
  GITHUB_ACTION=__run.
  GITHUB_ACTION_PATH=.
  GITHUB_ACTION_REPOSITORY=.
  GITHUB_ACTIONS=true.
  GITHUB_ACTOR=eric2003.
  GITHUB_ACTOR_ID=7891909.
  GITHUB_API_URL=https://api.github.com.
  GITHUB_BASE_REF=.
  GITHUB_ENV=/home/runner/work/_temp/_runner_file_commands/set_env_0b71d8bf-192c-41be-b943-63d02078de36.
  GITHUB_EVENT_NAME=push.
  GITHUB_EVENT_PATH=/home/runner/work/_temp/_github_workflow/event.json.
  GITHUB_GRAPHQL_URL=https://api.github.com/graphql.
  GITHUB_HEAD_REF=.
  GITHUB_JOB=Display-GitHub-Default-environment-variables.
  GITHUB_OUTPUT=/home/runner/work/_temp/_runner_file_commands/set_output_0b71d8bf-192c-41be-b943-63d02078de36.
  GITHUB_PATH=/home/runner/work/_temp/_runner_file_commands/add_path_0b71d8bf-192c-41be-b943-63d02078de36.
  GITHUB_REF=refs/heads/develop.
  GITHUB_REF_NAME=develop.
  GITHUB_REF_PROTECTED=false.
  GITHUB_REF_TYPE=branch.
  GITHUB_REPOSITORY=eric2003/ModernPrj.
  GITHUB_REPOSITORY_ID=735798203.
  GITHUB_REPOSITORY_OWNER=eric2003.
  GITHUB_REPOSITORY_OWNER_ID=7891909.
  GITHUB_RETENTION_DAYS=90.
  GITHUB_RUN_ATTEMPT=1.
  GITHUB_RUN_ID=7354347118.
  GITHUB_RUN_NUMBER=4.
  GITHUB_SERVER_URL=https://github.com.
  GITHUB_SHA=f62ec4d27ca3c57169a8f108fb46574cfa785098.
  GITHUB_STEP_SUMMARY=/home/runner/work/_temp/_runner_file_commands/step_summary_0b71d8bf-192c-41be-b943-63d02078de36.
  GITHUB_TRIGGERING_ACTOR=eric2003.
  GITHUB_WORKFLOW=GitHub Default environment variables.
  GITHUB_WORKFLOW_REF=eric2003/ModernPrj/.github/workflows/default-environment-variables.yml@refs/heads/develop.
  GITHUB_WORKFLOW_SHA=f62ec4d27ca3c57169a8f108fb46574cfa785098.
  GITHUB_WORKSPACE=/home/runner/work/ModernPrj/ModernPrj.
  RUNNER_ARCH=X64.
  RUNNER_DEBUG=.
  RUNNER_NAME=GitHub Actions 4.
  RUNNER_OS=Linux.
  RUNNER_TEMP=/home/runner/work/_temp.
  RUNNER_TOOL_CACHE=/opt/hostedtoolcache.
  
Github Project Multi Platform Info
::

  name: GitHub Environment Variables
  on: [push]
  jobs:
    my_job:
      runs-on: ${{matrix.os}}
      strategy:
        matrix:
          os: [ubuntu-latest, windows-latest]
      steps:
        - name: Default-environment-variables
          run: |    
            echo "github.action=${{ github.action }}."
            echo "github.action_path=${{ github.action_path }}."
            echo "github.action_ref=${{ github.action_ref }}."
            echo "github.action_repository=${{ github.action_repository }}."
            echo "github.action_status=${{ github.action_status }}."
            echo "github.actor=${{ github.actor }}."
            echo "github.actor_id=${{ github.actor_id }}."
            echo "github.api_url=${{ github.api_url }}."
            echo "github.base_ref=${{ github.base_ref }}."
            echo "github.env=${{ github.env }}."
            echo "github.event=${{ github.event }}."
            echo "github.event_name=${{ github.event_name }}."
            echo "github.event_path=${{ github.event_path }}."
            echo "github.graphql_url=${{ github.graphql_url }}."
            echo "github.head_ref=${{ github.head_ref }}."
            echo "github.job=${{ github.job }}."
            echo "github.path=${{ github.path }}."
            echo "github.ref=${{ github.ref }}."
            echo "github.ref_name=${{ github.ref_name }}."
            echo "github.ref_protected=${{ github.ref_protected }}."
            echo "github.ref_type=${{ github.ref_type }}."
            echo "github.repository=${{ github.repository }}."
            echo "github.repository_id=${{ github.repository_id }}."
            echo "github.repository_owner=${{ github.repository_owner }}."
            echo "github.repositoryUrl=${{ github.repositoryUrl }}."
            echo "github.retention_days=${{ github.retention_days }}."
            echo "github.run_id=${{ github.run_id }}."
            echo "github.run_number=${{ github.run_number }}."
            echo "github.run_attempt=${{ github.run_attempt }}."
            echo "github.secret_source=${{ github.secret_source }}."
            echo "github.server_url=${{ github.server_url }}."
            echo "github.sha=${{ github.sha }}."
            echo "github.token=${{ github.token }}."
            echo "github.triggering_actor=${{ github.triggering_actor }}."
            echo "github.workflow=${{ github.workflow }}."
            echo "github.workflow_ref=${{ github.workflow_ref }}."
            echo "github.workflow_sha=${{ github.workflow_sha }}."
            echo "github.workspace=${{ github.workspace }}."

Ubuntu Results:
::

  github.action=__run.
  github.action_path=.
  github.action_ref=.
  github.action_repository=.
  github.action_status=.
  github.actor=eric2003.
  github.actor_id=7891909.
  github.api_url=https://api.github.com.
  github.base_ref=.
  github.env=/home/runner/work/_temp/_runner_file_commands/set_env_20bde834-6a27-4b10-b3f6-789728256309.
  github.event=Object.
  github.event_name=push.
  github.event_path=/home/runner/work/_temp/_github_workflow/event.json.
  github.graphql_url=https://api.github.com/graphql.
  github.head_ref=.
  github.job=my_job.
  github.path=/home/runner/work/_temp/_runner_file_commands/add_path_20bde834-6a27-4b10-b3f6-789728256309.
  github.ref=refs/heads/develop.
  github.ref_name=develop.
  github.ref_protected=false.
  github.ref_type=branch.
  github.repository=eric2003/ModernPrj.
  github.repository_id=735798203.
  github.repository_owner=eric2003.
  github.repositoryUrl=git://github.com/eric2003/ModernPrj.git.
  github.retention_days=90.
  github.run_id=7355184823.
  github.run_number=9.
  github.run_attempt=1.
  github.secret_source=Actions.
  github.server_url=https://github.com.
  github.sha=99107ae8911cfc72716aecaf56232ad44606bae9.
  github.token=***.
  github.triggering_actor=eric2003.
  github.workflow=GitHub Environment Variables.
  github.workflow_ref=eric2003/ModernPrj/.github/workflows/default-environment-variables.yml@refs/heads/develop.
  github.workflow_sha=99107ae8911cfc72716aecaf56232ad44606bae9.
  github.workspace=/home/runner/work/ModernPrj/ModernPrj.

Windows Results:
::

  github.action=__run.
  github.action_path=.
  github.action_ref=.
  github.action_repository=.
  github.action_status=.
  github.actor=eric2003.
  github.actor_id=7891909.
  github.api_url=https://api.github.com.
  github.base_ref=.
  github.env=D:\a\_temp\_runner_file_commands\set_env_c7524d19-ed50-4919-b0ae-70648f47e68c.
  github.event=Object.
  github.event_name=push.
  github.event_path=D:\a\_temp\_github_workflow\event.json.
  github.graphql_url=https://api.github.com/graphql.
  github.head_ref=.
  github.job=my_job.
  github.path=D:\a\_temp\_runner_file_commands\add_path_c7524d19-ed50-4919-b0ae-70648f47e68c.
  github.ref=refs/heads/develop.
  github.ref_name=develop.
  github.ref_protected=false.
  github.ref_type=branch.
  github.repository=eric2003/ModernPrj.
  github.repository_id=735798203.
  github.repository_owner=eric2003.
  github.repositoryUrl=git://github.com/eric2003/ModernPrj.git.
  github.retention_days=90.
  github.run_id=7355184823.
  github.run_number=9.
  github.run_attempt=1.
  github.secret_source=Actions.
  github.server_url=https://github.com.
  github.sha=99107ae8911cfc72716aecaf56232ad44606bae9.
  github.token=***.
  github.triggering_actor=eric2003.
  github.workflow=GitHub Environment Variables.
  github.workflow_ref=eric2003/ModernPrj/.github/workflows/default-environment-variables.yml@refs/heads/develop.
  github.workflow_sha=99107ae8911cfc72716aecaf56232ad44606bae9.
  github.workspace=D:\a\ModernPrj\ModernPrj.

  
Run echo "$GITHUB_CONTEXT"
::

  {
    "token": "***",
    "job": "dump_contexts_to_log",
    "ref": "refs/heads/develop",
    "sha": "8d54f5ef375c67b2b7fee9a49f03213b96eba3f9",
    "repository": "eric2003/ModernPrj",
    "repository_owner": "eric2003",
    "repository_owner_id": "7891909",
    "repositoryUrl": "git://github.com/eric2003/ModernPrj.git",
    "run_id": "7357250199",
    "run_number": "1",
    "retention_days": "90",
    "run_attempt": "1",
    "artifact_cache_size_limit": "10",
    "repository_visibility": "public",
    "repo-self-hosted-runners-disabled": false,
    "enterprise-managed-business-id": "",
    "repository_id": "735798203",
    "actor_id": "7891909",
    "actor": "eric2003",
    "triggering_actor": "eric2003",
    "workflow": "Context testing",
    "head_ref": "",
    "base_ref": "",
    "event_name": "push",
    "event": {
      "after": "8d54f5ef375c67b2b7fee9a49f03213b96eba3f9",
      "base_ref": null,
      "before": "99107ae8911cfc72716aecaf56232ad44606bae9",
      "commits": [
        {
          "author": {
            "email": "fantasy_2003_@hotmail.com",
            "name": "eric",
            "username": "eric2003"
          },
          "committer": {
            "email": "fantasy_2003_@hotmail.com",
            "name": "eric",
            "username": "eric2003"
          },
          "distinct": true,
          "id": "8d54f5ef375c67b2b7fee9a49f03213b96eba3f9",
          "message": "add context.ymal",
          "timestamp": "2023-12-29T22:00:21+08:00",
          "tree_id": "4a07e48baec9177de003b537b74e3930c0416c65",
          "url": "https://github.com/eric2003/ModernPrj/commit/8d54f5ef375c67b2b7fee9a49f03213b96eba3f9"
        }
      ],
      "compare": "https://github.com/eric2003/ModernPrj/compare/99107ae8911c...8d54f5ef375c",
      "created": false,
      "deleted": false,
      "forced": false,
      "head_commit": {
        "author": {
          "email": "fantasy_2003_@hotmail.com",
          "name": "eric",
          "username": "eric2003"
        },
        "committer": {
          "email": "fantasy_2003_@hotmail.com",
          "name": "eric",
          "username": "eric2003"
        },
        "distinct": true,
        "id": "8d54f5ef375c67b2b7fee9a49f03213b96eba3f9",
        "message": "add context.ymal",
        "timestamp": "2023-12-29T22:00:21+08:00",
        "tree_id": "4a07e48baec9177de003b537b74e3930c0416c65",
        "url": "https://github.com/eric2003/ModernPrj/commit/8d54f5ef375c67b2b7fee9a49f03213b96eba3f9"
      },
      "pusher": {
        "email": "fantasy_2003_@hotmail.com",
        "name": "eric2003"
      },
      "ref": "refs/heads/develop",
      "repository": {
        "allow_forking": true,
        "archive_url": "https://api.github.com/repos/eric2003/ModernPrj/{archive_format}{/ref}",
        "archived": false,
        "assignees_url": "https://api.github.com/repos/eric2003/ModernPrj/assignees{/user}",
        "blobs_url": "https://api.github.com/repos/eric2003/ModernPrj/git/blobs{/sha}",
        "branches_url": "https://api.github.com/repos/eric2003/ModernPrj/branches{/branch}",
        "clone_url": "https://github.com/eric2003/ModernPrj.git",
        "collaborators_url": "https://api.github.com/repos/eric2003/ModernPrj/collaborators{/collaborator}",
        "comments_url": "https://api.github.com/repos/eric2003/ModernPrj/comments{/number}",
        "commits_url": "https://api.github.com/repos/eric2003/ModernPrj/commits{/sha}",
        "compare_url": "https://api.github.com/repos/eric2003/ModernPrj/compare/{base}...{head}",
        "contents_url": "https://api.github.com/repos/eric2003/ModernPrj/contents/{+path}",
        "contributors_url": "https://api.github.com/repos/eric2003/ModernPrj/contributors",
        "created_at": 1703569045,
        "default_branch": "main",
        "deployments_url": "https://api.github.com/repos/eric2003/ModernPrj/deployments",
        "description": null,
        "disabled": false,
        "downloads_url": "https://api.github.com/repos/eric2003/ModernPrj/downloads",
        "events_url": "https://api.github.com/repos/eric2003/ModernPrj/events",
        "fork": false,
        "forks": 0,
        "forks_count": 0,
        "forks_url": "https://api.github.com/repos/eric2003/ModernPrj/forks",
        "full_name": "eric2003/ModernPrj",
        "git_commits_url": "https://api.github.com/repos/eric2003/ModernPrj/git/commits{/sha}",
        "git_refs_url": "https://api.github.com/repos/eric2003/ModernPrj/git/refs{/sha}",
        "git_tags_url": "https://api.github.com/repos/eric2003/ModernPrj/git/tags{/sha}",
        "git_url": "git://github.com/eric2003/ModernPrj.git",
        "has_discussions": false,
        "has_downloads": true,
        "has_issues": true,
        "has_pages": false,
        "has_projects": true,
        "has_wiki": true,
        "homepage": null,
        "hooks_url": "https://api.github.com/repos/eric2003/ModernPrj/hooks",
        "html_url": "https://github.com/eric2003/ModernPrj",
        "id": 735798203,
        "is_template": false,
        "issue_comment_url": "https://api.github.com/repos/eric2003/ModernPrj/issues/comments{/number}",
        "issue_events_url": "https://api.github.com/repos/eric2003/ModernPrj/issues/events{/number}",
        "issues_url": "https://api.github.com/repos/eric2003/ModernPrj/issues{/number}",
        "keys_url": "https://api.github.com/repos/eric2003/ModernPrj/keys{/key_id}",
        "labels_url": "https://api.github.com/repos/eric2003/ModernPrj/labels{/name}",
        "language": "C++",
        "languages_url": "https://api.github.com/repos/eric2003/ModernPrj/languages",
        "license": {
          "key": "mit",
          "name": "MIT License",
          "node_id": "MDc6TGljZW5zZTEz",
          "spdx_id": "MIT",
          "url": "https://api.github.com/licenses/mit"
        },
        "master_branch": "main",
        "merges_url": "https://api.github.com/repos/eric2003/ModernPrj/merges",
        "milestones_url": "https://api.github.com/repos/eric2003/ModernPrj/milestones{/number}",
        "mirror_url": null,
        "name": "ModernPrj",
        "node_id": "R_kgDOK9tjuw",
        "notifications_url": "https://api.github.com/repos/eric2003/ModernPrj/notifications{?since,all,participating}",
        "open_issues": 0,
        "open_issues_count": 0,
        "owner": {
          "avatar_url": "https://avatars.githubusercontent.com/u/7891909?v=4",
          "email": "fantasy_2003_@hotmail.com",
          "events_url": "https://api.github.com/users/eric2003/events{/privacy}",
          "followers_url": "https://api.github.com/users/eric2003/followers",
          "following_url": "https://api.github.com/users/eric2003/following{/other_user}",
          "gists_url": "https://api.github.com/users/eric2003/gists{/gist_id}",
          "gravatar_id": "",
          "html_url": "https://github.com/eric2003",
          "id": 7891909,
          "login": "eric2003",
          "name": "eric2003",
          "node_id": "MDQ6VXNlcjc4OTE5MDk=",
          "organizations_url": "https://api.github.com/users/eric2003/orgs",
          "received_events_url": "https://api.github.com/users/eric2003/received_events",
          "repos_url": "https://api.github.com/users/eric2003/repos",
          "site_admin": false,
          "starred_url": "https://api.github.com/users/eric2003/starred{/owner}{/repo}",
          "subscriptions_url": "https://api.github.com/users/eric2003/subscriptions",
          "type": "User",
          "url": "https://api.github.com/users/eric2003"
        },
        "private": false,
        "pulls_url": "https://api.github.com/repos/eric2003/ModernPrj/pulls{/number}",
        "pushed_at": 1703858432,
        "releases_url": "https://api.github.com/repos/eric2003/ModernPrj/releases{/id}",
        "size": 17,
        "ssh_url": "git@github.com:eric2003/ModernPrj.git",
        "stargazers": 0,
        "stargazers_count": 0,
        "stargazers_url": "https://api.github.com/repos/eric2003/ModernPrj/stargazers",
        "statuses_url": "https://api.github.com/repos/eric2003/ModernPrj/statuses/{sha}",
        "subscribers_url": "https://api.github.com/repos/eric2003/ModernPrj/subscribers",
        "subscription_url": "https://api.github.com/repos/eric2003/ModernPrj/subscription",
        "svn_url": "https://github.com/eric2003/ModernPrj",
        "tags_url": "https://api.github.com/repos/eric2003/ModernPrj/tags",
        "teams_url": "https://api.github.com/repos/eric2003/ModernPrj/teams",
        "topics": [],
        "trees_url": "https://api.github.com/repos/eric2003/ModernPrj/git/trees{/sha}",
        "updated_at": "2023-12-27T10:36:17Z",
        "url": "https://github.com/eric2003/ModernPrj",
        "visibility": "public",
        "watchers": 0,
        "watchers_count": 0,
        "web_commit_signoff_required": false
      },
      "sender": {
        "avatar_url": "https://avatars.githubusercontent.com/u/7891909?v=4",
        "events_url": "https://api.github.com/users/eric2003/events{/privacy}",
        "followers_url": "https://api.github.com/users/eric2003/followers",
        "following_url": "https://api.github.com/users/eric2003/following{/other_user}",
        "gists_url": "https://api.github.com/users/eric2003/gists{/gist_id}",
        "gravatar_id": "",
        "html_url": "https://github.com/eric2003",
        "id": 7891909,
        "login": "eric2003",
        "node_id": "MDQ6VXNlcjc4OTE5MDk=",
        "organizations_url": "https://api.github.com/users/eric2003/orgs",
        "received_events_url": "https://api.github.com/users/eric2003/received_events",
        "repos_url": "https://api.github.com/users/eric2003/repos",
        "site_admin": false,
        "starred_url": "https://api.github.com/users/eric2003/starred{/owner}{/repo}",
        "subscriptions_url": "https://api.github.com/users/eric2003/subscriptions",
        "type": "User",
        "url": "https://api.github.com/users/eric2003"
      }
    },
    "server_url": "https://github.com",
    "api_url": "https://api.github.com",
    "graphql_url": "https://api.github.com/graphql",
    "ref_name": "develop",
    "ref_protected": false,
    "ref_type": "branch",
    "secret_source": "Actions",
    "workflow_ref": "eric2003/ModernPrj/.github/workflows/context.yml@refs/heads/develop",
    "workflow_sha": "8d54f5ef375c67b2b7fee9a49f03213b96eba3f9",
    "workspace": "/home/runner/work/ModernPrj/ModernPrj",
    "action": "__run"
  }  
  
Run echo "$JOB_CONTEXT"
::

  {
    "status": "success"
  }
  
Run echo "$STEPS_CONTEXT"
::

  {} 
  
  
Run echo "$RUNNER_CONTEXT"
::

  {
    "os": "Linux",
    "arch": "X64",
    "name": "GitHub Actions 5",
    "environment": "github-hosted",
    "tool_cache": "/opt/hostedtoolcache",
    "temp": "/home/runner/work/_temp",
    "workspace": "/home/runner/work/ModernPrj"
  }
  

Run echo "$STRATEGY_CONTEXT"
::

  {
    "fail-fast": true,
    "job-index": 0,
    "job-total": 1,
    "max-parallel": 1
  }
  
  
Run echo "$MATRIX_CONTEXT"
::

  null  
  
How To Pass Data Between Workflows in GitHub Actions?
--------------------------------------------------------
#. `How To Pass Data Between Workflows in GitHub Actions? <https://codersee.com/how-to-pass-data-between-workflows-in-github-actions/>`_  


REST API
-----------
#. `GitHub Actions Secrets <https://docs.github.com/en/rest/actions/secrets?apiVersion=2022-11-28/>`_  


GitHub Actions Windows MSYS2 MS-MPI
------------------------------------
#. `GitHub Actions Windows MSYS2 MS-MPI <https://www.scivision.dev/github-actions-windows-mpi-msys/>`_  


