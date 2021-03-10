import time
import os
import random
import subprocess
import logging
from multiprocessing import Process, JoinableQueue, Lock

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
            logging.info('{} got {}'.format(os.getpid(), query_cmd))
        try:
            subprocess.call(query_cmd, shell=True)
        except Exception as e:
            with lock:
                logging.info("Failed to execute " + query_cmd)
        queue.task_done()

def invoke_producer_consumer(cmd_list, consumer_ctr):
    print('Parent process staring...')
    # Create the Queue object
    queue = JoinableQueue()
    # Create a lock object to synchronize resource access
    lock = Lock()
    consumers = []

    print('Setting up consumers.')
    # Create consumer processes
    for i in range(consumer_ctr):
        p = Process(target=consumer, args=(queue, lock))
        consumers.append(p)
    # Start the producers and consumer
    # The Python VM will launch new independent processes for each Process object
    print('Starting up consumer sub-processes.')
    for c in consumers:
        c.start()

    print('Putting jobs in queue.')
    # Create our producer processes by passing the producer function and it's arguments
    for cmd_str in cmd_list:
        queue.put(cmd_str)

    # join() method that synchronizes our program
    queue.join()
    print('Stopping sub-processes and terminating.')
    for c in consumers:
        c.terminate()

    print('Parent process exiting...')



