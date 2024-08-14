#!/bin/bash

#docker login -u shenkers

is_dirty=false

if [ -n "$(git ls-files --others --exclude-standard)" ]; then
    echo "There are untracked files."
    is_dirty=true
fi

if ! git diff-index --quiet HEAD --; then
    echo "There are modified files"
    is_dirty=true
fi

if [ "$is_dirty" = true ]; then
    docker_tag=shenkers/cfdscore:dev
else
    docker_tag=shenkers/cfdscore
    git_commit_sha=$(git rev-parse --short HEAD)
fi

docker build -t ${docker_tag} ${git_commit_sha:+-t shenkers/cfdscore:$git_commit_sha} .

#docker push shenkers/cfdscore
