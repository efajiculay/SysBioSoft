#9grwg8eew3ew7rqw13525O35o54h56ayyeAeryseryteyhyrreydeyp
#pypi-AgENdGVzdC5weXBpLm9yZwIkMjhlODI0NjktOTNjOS00YjNjLWE5MzktNDE1NDZmYzA2NjdkAAIleyJwZXJtaXNzaW9ucyI6ICJ1c2VyIiwgInZlcnNpb24iOiAxfQAABiBNFUxK_KG7Hk3y0R55Xlu-NUjq3xqtSsyZ6p-f6dY1iQ

#https://test.pypi.org/project/BioSANS2020-efajiculay/0.0.1/
#pip install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay==0.0.1
#https://realpython.com/pypi-publish-python-package/


#https://realpython.com/pypi-publish-python-package/
python setup.py sdist bdist_wheel
#python setup.py --help-commands
#python setup.py sdist --help
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#twine upload dist/*


#https://realpython.com/python-import/#create-and-install-a-local-package
python -m pip install -e ./

#https://packaging.python.org/tutorials/packaging-projects/
#working
python -m build
python -m twine check dist/*
python -m twine upload --repository testpypi dist/* 

#working
python setup.py sdist bdist_wheel
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

#https://www.youtube.com/watch?v=9GfIDdIxhuc
#https://pyinstaller.readthedocs.io/en/stable/operating-mode.

#https://docs.anaconda.com/anaconda/install/silent-mode/
#https://godatadriven.com/blog/a-practical-guide-to-using-setup-py/

#pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple BioSANS2020==0.0.6
