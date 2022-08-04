from metabgc.src.producer_consumer import *
from metabgc.src.utils import *
import os
import subprocess
import metabgc.src.hmmrecord as hmmrecord
import re
import logging

"""
Function searches all FASTA file in a directory against a HMM in parallel. 
"""
def RunPCHMMDirectoryParallel(inputDir, hmmModel, sampleType, protType, window, interval, ouputDir, ncpus=4):
    logging.info('Number of pool processes:{0}.'.format(ncpus))
    hmmsearch_task_list = []
    for subdir, dirs, files in os.walk(inputDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if re.match(r".*.fasta$", file) and os.path.getsize(filePath) > 0:
                sampleStr = os.path.splitext(file)[0]
                hmm_task = hmmrecord.HMMTask(filePath,hmmModel,ouputDir,sampleType,sampleStr,protType, window, interval)
                hmmsearch_task_list.append(hmm_task)


    print('HMMER searching staring for: ' + hmmModel)
    # Create the Queue object
    queue = JoinableQueue()
    # Create a lock object to synchronize resource access
    lock = Lock()
    consumers = []

    print('Setting up consumers.')
    # Create consumer processes
    for i in range(ncpus):
        p = Process(target=HMMSearchConsumer, args=(queue, lock))
        consumers.append(p)
    # Start the producers and consumer
    # The Python VM will launch new independent processes for each Process object
    print('Starting up consumer sub-processes.')
    for c in consumers:
        c.start()

    print('Putting jobs in queue.')
    # Create our producer processes by passing the producer function and it's arguments
    for task_obj in hmmsearch_task_list:
        queue.put(task_obj)

    # join() method that synchronizes our program
    queue.join()
    print('Stopping sub-processes and terminating.')
    for c in consumers:
        c.terminate()
    print('HMMER searching exiting for: ' + hmmModel)

"""
Function searches FASTA file against HMM. 
"""
def HMMSearchConsumer(queue, lock):

    # Synchronize access to the console
    with lock:
        logging.info('Starting consumer => {}'.format(os.getpid()))
    # Run indefinitely
    while True:
        time.sleep(random.randint(0, 10))
        # If the queue is empty, queue.get() will block until the queue has data
        hmmsearch_task = queue.get()
        # Synchronize access to the console
        with lock:
            logging.info('{} got {}'.format(os.getpid(), hmmsearch_task.fastaFile))

        hmmTblFileName = hmmsearch_task.sampleStr + "__" + hmmsearch_task.interval + ".tbl"
        hmmTblFilePath = os.path.join(hmmsearch_task.ouputDir, hmmTblFileName)
        if not os.path.exists(hmmTblFilePath):
            cmd = "hmmsearch --cpu " + str(hmmsearch_task.ncpus) + " --F1 0.02 --F2 0.02 --F3 0.02 --tblout " + hmmTblFilePath + " " + hmmsearch_task.hmmFile + " "+ hmmsearch_task.fastaFile + " > /dev/null"
            with lock:
                logging.info('Running HMM Search with {0} against {1}.'.format(hmmsearch_task.fastaFile, hmmsearch_task.hmmFile))
                logging.info(cmd)
            subprocess.call(cmd, shell=True)
        else:
            with lock:
                logging.info("HMM search skipped... Using existing result for: " + hmmsearch_task.fastaFile)
        result_dict = parseHMM(hmmTblFilePath, "hmmer3-tab", hmmsearch_task.sampleType, hmmsearch_task.sampleStr,
                               hmmsearch_task.protType, hmmsearch_task.window, hmmsearch_task.interval)
        hmmSearchFileName = hmmsearch_task.sampleStr + "__" + hmmsearch_task.interval + ".txt"
        hmmSearchFilePath = os.path.join(hmmsearch_task.ouputDir, hmmSearchFileName)
        createPandaDF(result_dict, hmmSearchFilePath)
        with lock:
            logging.info("Done Running HMM search with:" + hmmsearch_task.fastaFile)

        queue.task_done()
