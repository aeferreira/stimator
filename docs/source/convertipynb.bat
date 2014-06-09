
ipython nbconvert --to rst --template non_exec_rst.tpl solving.ipynb
ipython nbconvert --to rst --template executed_rst.tpl solving_executed.ipynb
copy solving_executed_files\*.* .\_static\solving

python fix_links.py solving_executed.rst
