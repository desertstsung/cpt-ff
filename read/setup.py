from distutils.core import setup, Extension

def main():
	setup(
		name="pycpt",
		version="0.1.0",
		description="Python interface for reading cpt file format file",
		author="Jay Tsung",
		author_email="dongjt@proton.me",
		ext_modules=[Extension("pycpt", ["readcpt_py.c"])])

if __name__ == "__main__":
	main()

