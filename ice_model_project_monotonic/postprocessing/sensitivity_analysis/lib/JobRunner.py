#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess		# import the module for calling external programs (creating subprocesses)
import sys
import os
import shutil
import threading,queue
import concurrent.futures


class JobRunner:
	"""
	This class encapsulates a thread queue, a runner function and provides a convenient
	interface for running jobs in parallel.
	"""
	def __init__(self, N):
		# create a queue for task scheduling
		# self.jobQueue = queue.Queue()
		self.N = N

	# here we define our solver thread function
	def process_function(self, job):
		"""This function runs a simulation solver."""
		try:
			subprocess.run(job, check=True, stdout=subprocess.PIPE)
			print(f"Process complete: {job}")
		except Exception as e:
			print(str(e))
		# # run as long as we have jobs in the Queue
		# while 1:
		# 	# get a new task from the Queue
		# 	job = self.jobQueue.get()
		# 	if job is None:
		# 		break
		# 	try:
		#		solverProcess = subprocess.Popen(job, stdout=subprocess.PIPE)
		# 		solverProcess.wait()
		# 	except Exception as e:
		# 		print(str(e))
		# 	# tell the Queue that the task was done
		# 	self.jobQueue.task_done()

	def run(self, jobs):
		with concurrent.futures.ProcessPoolExecutor(max_workers=self.N) as executor:
			executor.map(self.process_function, jobs)
		return
	# def run(self, jobs):
	# 	# now create as many threads as we need
	# 	for i in range(self.N):
	# 		t = threading.Thread(target=JobRunner.thread_function,args=[self])
	# 		t.daemon = True
	# 		t.start()
	# 	# now add project files to the queue
	# 	for job in jobs:
	# 		self.jobQueue.put(job)
	# 	# wait until all tasks are done
	# 	self.jobQueue.join()
	# 	return
