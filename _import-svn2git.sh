#!/bin/bash

# import the svn repo as git repo,
git svn clone -T trunk -b branches -t tags http://iveselyk@merlin.fit.vutbr.cz/svn/STK

# make the branches visible,
cd STK
git branch -a # show all branches, 
git checkout -b 64bit remotes/origin/64bit 
git checkout -b STree-devel remotes/origin/STree-devel
git checkout -b latgen remotes/origin/latgen
git checkout -b rd_xforms remotes/origin/rd_xforms
git checkout -b tags_release-2.0_alpha1 remotes/origin/tags/release-2.0_alpha1
git checkout -b tags_v0.0.0 remotes/origin/tags/v0.0.0

# add the location of remote repo,
git remote add origin https://vesis84@github.com/stk-developers/stk.git

# push the branches into the remote repository,
for branch in master 64bit STree-devel latgen rd_xforms tags_release-2.0_alpha1 tags_v0.0.0; do
  git checkout $branch
  git push -u origin $branch
done

