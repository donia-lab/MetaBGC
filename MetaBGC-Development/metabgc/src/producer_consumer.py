import time
import os
import random
import subprocess
from multiprocessing import Process, Queue, Lock

# Producer function that places data on the Queue
def producer(queue, lock, jobs):
    # Synchronize access to the console
    with lock:
        print('Starting producer => {}'.format(os.getpid()))
    # Place our names on the Queue
    for job_str in jobs:
        queue.put(job_str)
    # Synchronize access to the console
    with lock:
        print('Producer {} exiting...'.format(os.getpid()))

# The consumer function takes data off of the Queue and runs the command line
def consumer(queue, lock):
    # Synchronize access to the console
    with lock:
        print('Starting consumer => {}'.format(os.getpid()))
    # Run indefinitely
    while True:
        time.sleep(random.randint(0, 10))
        # If the queue is empty, queue.get() will block until the queue has data
        query_cmd = queue.get()
        # Synchronize access to the console
        with lock:
            print('{} got {}'.format(os.getpid(), query_cmd))
        subprocess.call(query_cmd, shell=True)

def invoke_producer_consumer(cmd_list, consumer_ctr):
    print('Parent process staring...')
    # Create the Queue object
    queue = Queue()
    # Create a lock object to synchronize resource access
    lock = Lock()
    queue = Queue()
    producers = []
    consumers = []

    print('Setting up producers.')
    # Create our producer processes by passing the producer function and it's arguments
    producers.append(Process(target=producer, args=(queue, lock, cmd_list)))

    print('Setting up consumers.')
    # Create consumer processes
    for i in range(consumer_ctr):
        p = Process(target=consumer, args=(queue, lock))
        # This is critical! The consumer function has an infinite loop
        # Which means it will never exit unless we set daemon to true
        p.daemon = True
        consumers.append(p)
    # Start the producers and consumer
    # The Python VM will launch new independent processes for each Process object

    print('Starting up sub-processes.')
    for p in producers:
        p.start()
    for c in consumers:
        c.start()

    # Like threading, we have a join() method that synchronizes our program
    for p in producers:
        p.join()

    print('Parent process exiting...')



