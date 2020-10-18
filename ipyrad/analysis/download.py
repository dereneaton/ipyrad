#!/usr/bin/env python


import os
import gzip
import requests



class Download:
	"""
	Download a large file by streaming in chunks using the requests library.
	The file from 'url' is downloaded to the location 'path'. The file will
	not be downloaded if it already exists at this path, unless you overwrite
	using the force option.
	"""
	def __init__(self, url, path, gunzip=False, force=False):
		self.url = url
		self.path = path
		self.force = force
		self.gunzip = gunzip

		self.run()
		if self.gunzip:
			try:
				self.gunzip_file()
			except Exception:
				print("error: couldn't gunzip file.")


	def run(self):
		# only run if the reference doesn't already exist
		if (not os.path.exists(self.path)) and (not self.force):
	
			# open a stream to url and write to file 1Mb at a time.
			res = requests.get(self.url, stream=True)
			with open(self.path, 'wb') as out:
				for chunk in res.iter_content(chunk_size=1024*1024):
					if chunk:
						out.write(chunk)
			print("successful download: {}".format(self.path))
		else:
			print("file already exists")


	def gunzip_file(self):
		# make a decompressed copy of the file
		rname = self.path.split(".gz")[0]
		if not os.path.exists(rname):
			with open(rname, 'w') as out:
				out.write(gzip.open(self.path).read().decode())