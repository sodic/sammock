# pre-commit.sh
STASH_NAME="pre-commit-$(date +%s)"
git stash save -q --keep-index $STASH_NAME

# Test prospective commit
echo 'Running tests...'
python3 -m unittest test_sammock.py
if [ $? -ne 0 ]; then
    echo "All tests need to pass in order to commit."
    exit 1
fi

STASHES=$(git stash list)
if [[ $STASHES == "$STASH_NAME" ]]; then
  git stash pop -q
fi
