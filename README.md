# scripts
SBB - Desktop

# initially!
cd /Users/user/Desktop/GitHub_jshinb/scripts

git add 2015-10-28_draw_pen_likelihood_function.R
git commit -m "drawing likelihood functions for the standard and penalized logistic regression - revisited"
git status
git log

# after create a repository on the github website! [run all three lines again??]
git commit -m "Change READNE.md"

# run all these command lines again after committing.
git remote add origin https://github.com/jshinb/scripts.git
git remote -v
git push -u origin master
# got the error messages!

--------------------------------------------------
To https://github.com/jshinb/scripts.git
 ! [rejected]        master -> master (non-fast-forward)
error: failed to push some refs to 'https://github.com/jshinb/scripts.git'
hint: Updates were rejected because the tip of your current branch is behind
hint: its remote counterpart. Merge the remote changes (e.g. 'git pull')
hint: before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
--------------------------------------------------
