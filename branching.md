# Below is a brief guide to branching.

## Create a new branch:

git checkout -b <branchname>
Always create branches from the latest version of main

# Switch branches:

git checkout <branchname>

# List all branches:

git branch

# branchname conventions: 

describe the feature contributed in the branch

# merging branches:

git push as usual. Go to the pull requests tab and create a request. If the branching strategy has been followed corectly, no merge conflicts should appear. If there are, manually resolve them. After creating the pull request, let another member briefly review the code before finalizing the merging.