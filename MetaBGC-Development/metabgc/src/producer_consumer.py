import time
import os
import random
from multiprocessing import Process, Queue, Lock


# Producer function that places data on the Queue
def producer(queue, lock, names):
    # Synchronize access to the console
    with lock:
        print('Starting producer => {}'.format(os.getpid()))

    # Place our names on the Queue
    for job_str in jobs:
        time.sleep(random.randint(0, 10))
        queue.put(job)

    # Synchronize access to the console
    with lock:
        print('Producer {} exiting...'.format(os.getpid()))


# The consumer function takes data off of the Queue
def consumer(queue, lock):
    # Synchronize access to the console
    with lock:
        print('Starting consumer => {}'.format(os.getpid()))

    # Run indefinitely
    while True:
        time.sleep(random.randint(0, 10))

        # If the queue is empty, queue.get() will block until the queue has data
        name = queue.get()

        # Synchronize access to the console
        with lock:
            print('{} got {}'.format(os.getpid(), name))

