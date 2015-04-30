import random
import threading


def list_append(count, id, out_list):
    """
    Creates an empty list and then appends a 
    random number to the list 'count' number
    of times. A CPU-heavy operation!
    """
    for i in range(count):
        out_list.append(random.random())

if __name__ == "__main__":
    size = 10000000   # Number of random numbers to add
    threads = 4   # Number of threads to create
    
    # Create a list of jobs and then iterate through
    # the number of threads appending each thread to
    # the job list 
    jobs = []
    for i in range(0, threads):
        out_list = list()
        thread = threading.Thread(target=list_append(size, i, out_list))
        jobs.append(thread)

    # Start the threads (i.e. calculate the random number lists)
    for j in jobs:
        j.start()

    # Ensure all of the threads have finished
    for j in jobs:
        j.join()

    print "List processing complete."
