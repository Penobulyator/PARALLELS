#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

int myTasksNumber;
int *myTasks;

int othersTaskNumber;
int *othersTasksQueue;
int othersTaskTail = 0;
int othersTaskHead = 0;

int leave = 0;

#define TASK_DONE -1

int rank;
int size;

pthread_mutex_t mutex;

#define MPI_RESPONSE_TAG 1
int noTasksResponse = -12345; //any number to send as a "no tasks response"


void doMyTasks() {
	for (int i = 0; i < myTasksNumber; i++) {
		//lock mutex so SEND thread doesn't access tasks
		pthread_mutex_lock(&mutex);
		if (myTasks[i] == TASK_DONE) {
			//SEND thread did this task for us, so we will just skip it
			printf("Process number %d did task %d / %d (some other process did it for me)\n", rank, i + 1, myTasksNumber);
			pthread_mutex_unlock(&mutex);
		}
		else {
			int time = myTasks[i];
			myTasks[i] = TASK_DONE;

			pthread_mutex_unlock(&mutex);

			sleep(time);
			printf("Process number %d did task %d / %d\n", rank, i + 1, myTasksNumber);
		}
	}
}
void doOthersTasks() {
	while (1) {
		pthread_mutex_lock(&mutex);
		if (othersTaskTail != othersTaskHead) {
			//if we have some tasks in others tasks queue
			int time = othersTasksQueue[othersTaskTail];
			othersTaskTail++;

			pthread_mutex_unlock(&mutex);
			sleep(time);
			printf("Process number %d did task of weigth %d for some other process\n", rank, time);
		}
		else {
			if (leave == 1)
			{
				pthread_mutex_unlock(&mutex);
				return;
			}
			else {
				pthread_mutex_unlock(&mutex);
				sleep(2);
			}
		}
	}
}
void* doThread(void* args){
	doMyTasks();

	doOthersTasks();
	//try to get some more tasks from other processes
	
	//syncronize with other process
	MPI_Barrier(MPI_COMM_WORLD);

	//tell SEND thread that we are done
	printf("Process number %d done his job\n", rank);

	return NULL;
}

int countLoad() {
	//pthread_mutex_lock(&mutex);
	int weigthSum = 0;
	for (int i = 0; i < myTasksNumber; i++) {
		if (myTasks[i] != TASK_DONE)
			weigthSum += myTasks[i];
	}

	for (int i = othersTaskTail; i < othersTaskHead; i++) {
		if (othersTasksQueue[i] != TASK_DONE)
			weigthSum += othersTasksQueue[i];
	}

	//pthread_mutex_unlock(&mutex);
	return weigthSum;
}

void receiveTask(int processNumber) {
	pthread_mutex_lock(&mutex);

	int answer = 0;
	MPI_Recv(&answer, 1, MPI_INT, processNumber, MPI_RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	if (answer != noTasksResponse) {
		//printf("Process number %d received task from process number %d\n", rank);
		othersTasksQueue[othersTaskHead] = answer;
		othersTaskHead++;
		//printf("Process number %d received task\n", rank);
	}
	else {
		//printf("Process number %d received NO_TASKS_RESPONSE\n", rank);
	}

	pthread_mutex_unlock(&mutex);
	
}

void sendTask(int processNumber) {

	//send task to process with lowest load
	pthread_mutex_lock(&mutex); 

	int taskSendedFlag = 0;
	for (int i = 0; i < myTasksNumber; i++)
		//find some undone task
		if (myTasks[i] != TASK_DONE) {
			int answer = myTasks[i];
			myTasks[i] = TASK_DONE;

			MPI_Send(&answer, 1, MPI_INT, processNumber, MPI_RESPONSE_TAG, MPI_COMM_WORLD);
			pthread_mutex_unlock(&mutex);

			//printf("Process number %d sending task to process %d\n", rank, processNumber);
			taskSendedFlag = 1;
			break;
		}
	//if no task to send was found
	if (taskSendedFlag == 0)
	{
		//printf("Process number %d sending NO_TASKS_RESPONSE\n", rank);
		MPI_Send(&noTasksResponse, 1, MPI_INT, processNumber, MPI_RESPONSE_TAG, MPI_COMM_WORLD);
		pthread_mutex_unlock(&mutex);
	}
}

void* sendThread(void* args) {
	while (1) {
		sleep(10);
		MPI_Barrier(MPI_COMM_WORLD);

		//find process with lowest load
		int load[2];
		load[0] = countLoad();
		load[1] = rank;

		MPI_Barrier(MPI_COMM_WORLD);
		printf("\t Current load for process number %d: %d\n", rank, load[0]);

		int minLoad[2];
		int maxLoad[2];
		MPI_Allreduce(load, &minLoad, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
		MPI_Allreduce(	load, &maxLoad, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

		if (minLoad[0] == maxLoad[0])
		{
			pthread_mutex_lock(&mutex);
			leave = 1;
			pthread_mutex_unlock(&mutex);
			printf("Process number %d ending SEND thread\n", rank);
			return NULL;
		}

		if (minLoad[1] == rank) {
			//if we are process with lowest load
			printf("Process number %d has lowes load\n", minLoad[1]);
			receiveTask(maxLoad[1]);
		}
		else if (maxLoad[1] == rank) {
			//if we are process with highest
			printf("Process number %d has highest load\n", maxLoad[1]);
			sendTask(minLoad[1]);
		}
	}
}

int createThreads() {
	//init threads attributes
	pthread_attr_t attrs;
	if (pthread_attr_init(&attrs) != 0) {
		perror("Cannot initialize attributes!");
		return 1;
	};

	//set threads joinable
	if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE) != 0) {
		perror("Error in setting attributes!");
		return 1;
	}

	//create DO, SEND threads
	pthread_t thrs[2];
	if (pthread_create(&thrs[0], &attrs, &doThread, NULL) != 0) {
		perror("Cannot create a DO thread!");
		return 1;
	}
	if (pthread_create(&thrs[1], &attrs, &sendThread, NULL) != 0) {
		perror("Cannot create a SEND thread!");
		return 1;
	}

	//free attributes
	pthread_attr_destroy(&attrs);

	//wait for DO and SEND threads to end threir job
	if (pthread_join(thrs[0], NULL) != 0) {
		perror("Cannot join a thread");
		return 1;
	}
	if (pthread_join(thrs[1], NULL) != 0) {
		perror("Cannot join a thread");
		return 1;
	}

	return 0;
}

void createMyTaskList() {
	myTasksNumber = 12 / (rank + 1);
	myTasks = (int*)malloc(myTasksNumber * sizeof(myTasks));
	for (int i = 0; i < myTasksNumber; i++)
		myTasks[i] = 5 + abs(rank - (i%size));

}

void createOthersTaskQueue() {
	othersTaskNumber = 10000;
	othersTasksQueue = (int*)calloc(othersTaskNumber, sizeof(othersTaskNumber));

	othersTaskTail = 0;
	othersTaskHead = 0;
}

int main(int argc, char** argv) {
	//init thread usage
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided != MPI_THREAD_MULTIPLE) {
		printf("%d \n", provided);
		MPI_Finalize();
		return 1;
	}

	//get rank and size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//init mutex
	pthread_mutex_init(&mutex, NULL);

	//init tasks
	createMyTaskList();
	createOthersTaskQueue();

	createThreads();

	printf("Programm done for process number %d\n", rank);

	return 0;

}