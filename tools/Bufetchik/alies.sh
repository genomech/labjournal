alias me-commit='msg="["$(date -R)"] "$(whoami)" commit at "$(hostnamectl --static)""; git add *; git commit -a -m "$msg"; git push'
