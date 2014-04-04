//
//  main.c
//  2048
//
//  Created by Timothy Johnson on 4/1/14.
//  Copyright (c) 2014 Timothy Johnson. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#define _size 4
#define _minProb 0.0001
#define _maxMoves 1000000
#define _possMoves 4
#define _numTests 10
#define _timeLimit 100000.0
#define _maxRowPoss 65536
#define _maxRank 16
#define _maxDepth 8

void expMove(unsigned long board, float expScore[], float currProb, int currDepth);
float expMove_low(unsigned long board, float currProb, int currDepth);
unsigned long doMove(unsigned long board, char move);
static inline int numEmpty(unsigned long board);
float max(float array[], int size);
int intLog(int val, int base);
int maxIndex_float(float array[], int size);
int maxIndex_int(int array[], int size);
unsigned long fillSquare(unsigned long board, int fillIndex, int toFill);
unsigned long updateBoard(unsigned long board, char move);
void printBoard_long(unsigned long board);
void printBoard_array(int board[][_size]);
void runTest();
char nextMove(unsigned long board);
char newNextMove(int powerBoard[][_size]);
float heuristic(unsigned long board);
static inline unsigned long transpose(unsigned long board);
void init_tables();
float calcRowHeur(int line[]);
unsigned long newLogBoard(int board[][_size]);
void newPowerBoard(unsigned long board, int newBoard[][_size]);
static inline int scoreVal(int boardVal);
int scoreBoard(unsigned long board);
int reverseRow(int row);

int row_lshift_table[_maxRowPoss];
int row_rshift_table[_maxRowPoss];
int row_score_table[_maxRowPoss];
int row_heur_table[_maxRowPoss];
int intTwoPower[_maxRank];
int offset[_size][_size];
int max_table[_maxRowPoss];
int numEmpty_table[_maxRowPoss];
long transpose_table[4][_maxRowPoss];
char allMoves[_possMoves] = {'U','D','L','R'};

int fillOptions = 2;
int toFill[2] = {1, 2}; //Since we are using logs, we are filling with 2 = 2**1, and 4 = 2**2
float prob[2] = {0.75, 0.25};
int expanded;
int currScore;

int main()
{
	clock_t t5, t6;
	t5 = clock();
	init_tables();
	t6 = clock();
	printf("Init time: %f\n",((float)t6 - (float)t5)/1000000);
	//runTest();
	for(int rep = 1; rep <= _numTests; rep++)
	{
		clock_t t1, t2, t3, t4;
		t1 = clock();
		printf("Starting: %d\n", rep);
		srand((int)time(NULL) + rep);
        
		unsigned long board = 0;
		unsigned long newBoard = 0;
		bool cont = true;
		int currMoveNum = 0;
		float maxMoveTime = -1;
		while(cont && currMoveNum < _maxMoves && maxMoveTime < _timeLimit)
		{
			t3 = clock();
			
			currMoveNum++;
            expanded = 0;
            /*int powerBoard[_size][_size];
            newPowerBoard(board, powerBoard);
			char move1 = newNextMove(powerBoard);
            char move2 = nextMove(board);
            if(move1 == move2)
                newBoard = updateBoard(board, move1);
            else printf("Error! Moves  do not match.\n");*/
            char move1 = nextMove(board);
            newBoard = updateBoard(board, move1);
			if(newBoard != board)
			{
				cont = true;
				board = newBoard;
                currScore = scoreBoard(board);
				
				t4 = clock();
				float timeElapsed = ((float)t4 - (float)t3)/1000;
				if(timeElapsed > maxMoveTime)
					maxMoveTime = timeElapsed;
				
                /*printf("Nodes expanded: %d\n",expanded);
				printf("Move %d: %c, %f ms\n",currMoveNum,move1,timeElapsed);
                printf("Time per node expanded: %f ns\n",timeElapsed*1000000/expanded);
                printf("Heuristic score: %f\n", heuristic(board));
				printBoard_long(board);*/
			}
			else
			{
				cont = false;
				//printf("Game over!\n");
			}
		}
		t2 = clock();
		printf("CPU time elapsed: %f s\n", ((float)t2 - (float)t1)/1000000);
		printf("Maximum move time: %f ms\n", maxMoveTime);
        printBoard_long(board);
	}
}

void init_tables()
{
    intTwoPower[0] = 0;
	int val = 2;
	for(int i = 1; i < _maxRank; i++)
	{
		intTwoPower[i] = val;
		val *= 2;
	}
    
	//Precompute a left shift
	int filled, locked;
	for(int i = 0; i < _maxRowPoss; i++)
	{
		int line[4] = {(i >> 12) & 0x0f, (i >> 8) & 0x0f, (i >> 4) & 0x0f, i & 0x0f};
		int newLine[4];
		memset(newLine, 0, 4*sizeof(int));
        
		filled = -1;
		locked = -1;
		for(int j = 0; j < 4; j++)
		{
			if(line[j] > 0 && line[j] == newLine[filled] && filled > locked)
			{
				newLine[filled]++;
				locked = filled;
			}
			else
			{
				newLine[filled + 1] = line[j];
				if(line[j] > 0)
					filled++;
			}
		}
		//printf("%d,%d,%d,%d\n",line[0],line[1],line[2],line[3]);
		//printf("%d,%d,%d,%d\n\n",newLine[0],newLine[1],newLine[2],newLine[3]);
        int newVal = (newLine[0] << 12) + (newLine[1] << 8) + (newLine[2] << 4) + newLine[3];
		row_lshift_table[i] = newVal;
        row_rshift_table[reverseRow(i)] = reverseRow(newVal);
		row_score_table[i] = scoreVal(line[0]) + scoreVal(line[1]) + scoreVal(line[2]) + scoreVal(line[3]);
        //printf("%d\n",row_score_table[i]);
		row_heur_table[i] = calcRowHeur(line);
        max_table[i] = maxIndex_int(line, _size);
        
        for(int j = 0; j < _size; j++)
        {
            transpose_table[j][i] = ((unsigned long)line[3] + ((unsigned long)line[2] << 16) + ((unsigned long)line[1] << 32) + ((unsigned long)line[0] << 48)) << (4*j);
        }
        
        for(int j = 0; j < _size; j++)
        {
            if(line[j] == 0)
                numEmpty_table[i]++;
        }
	}
	
	for(int i = 0; i < _size; i++)
	{
		for(int j = 0; j < _size; j++)
			offset[i][j] = sizeof(int)*(_size*(3 - i) + (3 - j));
	}
}

int reverseRow(int row)
{
    int newRow = 0;
    for(int i = 0; i < _size; i++)
    {
        int val = ((row >> 4*i) & 0x0f);
        newRow += val << 4*(3 - i);
    }
    return newRow;
}

void runTest()
{
	int board[_size][_size] = {{4, 8, 2, 2},{2, 2048, 64, 64},{256, 512, 128, 1024},{8, 2, 256, 0}};
	unsigned long logBoard = newLogBoard(board);
    printBoard_long(logBoard);
    printBoard_long(transpose(logBoard));
    //char move = nextMove(logBoard);
}

static inline int scoreVal(int boardVal)
{
	return (boardVal - 1)*intTwoPower[boardVal];
}

int scoreBoard(unsigned long board)
{
    int res = 0;
    int line;
    for(int i = 0; i < _size; i++)
    {
        line = (board >> (i << 4)) & 0xffff;
        res += row_score_table[line];
    }
    return res;
}

int intLog(int val, int base)
{
	if(val == 0)
		return 0;
	else
	{
		int res = 0;
		while(val > 1)
		{
			val /= base;
			res++;
		}
		return res;
	}
}

float calcRowHeur(int line[])
{
	float heurScore = 0.0;
	float heurEmpty = 100.0;
	float heurMaxEdge = 100.0;
	float heurNear = 25.0;
	float heurOrdered = 100.0;
	
	//Add for each empty square
	for(int i = 0; i < _size; i++)
	{
        if(line[i] == 0)
            heurScore += heurEmpty;
	}
	
	//Add for each maximum value in each row that is on the edge
	int maxCol = maxIndex_int(line, _size);
	if(maxCol == 0 || maxCol == _size - 1)
		heurScore += heurMaxEdge;
	
	//Add for each pair of adjacent values in each row that are exactly one apart
	for(int i = 0; i < _size - 1; i++)
	{
		if(line[i] == line[i + 1] + 1 || line[i] == line[i + 1] - 1)
			heurScore += heurNear;
	}
	
	//Add for each ordered row
	bool increasing = true;
	bool decreasing = true;
	for(int i = 0; i < _size - 1; i++)
	{
		if(line[i] <= line[i + 1])
			decreasing = false;
		if(line[i] >= line[i + 1])
			increasing = false;
	}
	if(increasing || decreasing)
		heurScore += heurOrdered;
    
	return heurScore;
}

unsigned long newLogBoard(int board[][_size])
{
	unsigned long newBoard = 0;
	unsigned long logVal;
	for(int i = 0; i < _size; i++)
	{
		for(int j = 0; j < _size; j++)
		{
			if(board[i][j] > 0)
			{
				logVal = intLog(board[i][j], 2);
				newBoard = newBoard | (logVal << offset[i][j]);
				//printf("%d,%d,%d,%ld,%ld\n",i,j,offset[i][j],logVal,newBoard);
			}
		}
	}
	return newBoard;
}

void newPowerBoard(unsigned long board, int newBoard[][_size])
{
	for(int i = 0; i < _size; i++)
	{
		for(int j = 0; j < _size; j++)
			newBoard[i][j] = intTwoPower[(board >> offset[i][j]) & 0x0f];
	}
}

char nextMove(unsigned long board)
{
	float expScore[_possMoves];
	memset(expScore, 0, _possMoves*sizeof(float));
	expMove(board, expScore, 1, 1);
	/*for(int i = 0; i < _possMoves; i++)
        printf("%f, ",expScore[i]);
    printf("\n");*/
	int bestMoveIndex = maxIndex_float(expScore, _possMoves);
	return allMoves[bestMoveIndex];
}

char newNextMove(int powerBoard[][_size])
{
  	init_tables();
  	long board = newLogBoard(powerBoard);
	char allMoves[_possMoves] = {'U','D','L','R'};
	float expScore[_possMoves];
	for(int i = 0; i < _possMoves; i++)
		expScore[i] = 0.0;
	expMove(board, expScore, 1, 1);
	
	int bestMoveIndex = maxIndex_float(expScore, _possMoves);
	char move = allMoves[bestMoveIndex];
    
  	/*if(move == 'U')
     printf("%s\n", "UP");
     else if(move == 'D')
     printf("%s\n", "DOWN");
     else if(move == 'L')
     printf("%s\n", "LEFT");
     else if(move == 'R')
     printf("%s\n", "RIGHT");*/
    return move;
}

static inline int numEmpty(unsigned long board)
{
    return numEmpty_table[board & 0xffff] + numEmpty_table[(board >> 16) & 0xffff]
            + numEmpty_table[(board >> 32) & 0xffff] + numEmpty_table[(board >> 48) & 0xffff];
}

float max(float array[], int size)
{
	float res = array[0];
	for(int i = 1; i < size; i++)
	{
		if(array[i] > res)
			res = array[i];
	}
	return res;
}

int maxIndex_float(float array[], int size)
{
	int res = 0;
	float best = array[0];
	for(int i = 1; i < size; i++)
	{
		if(array[i] > best)
		{
			best = array[i];
			res = i;
		}
	}
	return res;
}

int maxIndex_int(int array[], int size)
{
	int res = 0;
	int best = array[0];
	for(int i = 1; i < size; i++)
	{
		if(array[i] > best)
		{
			best = array[i];
			res = i;
		}
	}
	return res;
}

unsigned long fillSquare(unsigned long prevBoard, int fillIndex, int toFill)
{
	int count = 0;
	unsigned long newBoard = prevBoard;
	//printf("Filling: %ld\n",newBoard);
	for(int i = 0; i < _size && count <= fillIndex; i++)
	{
		for(int j = 0; j < _size; j++)
		{
			if(((prevBoard >> offset[i][j]) & 0x0f) == 0)
			{
				if(count == fillIndex)
					newBoard += ((unsigned long)toFill) << offset[i][j];
				
				count++;
			}
		}
	}
	//printf("Filled: %ld\n",newBoard);
	return newBoard;
}

unsigned long updateBoard(unsigned long board, char move)
{
	board = doMove(board, move);
	//printf("Updated board:\n");
	//printBoard_long(board);
	int empty = numEmpty(board);
	//printf("Empty: %d\n",empty);
	if(empty > 0)
	{
		int fillIndex = rand()%empty;
		int fillRand = rand()%4;
		//printf("fillRand: %d\n",fillRand);
		//We must pass the log of the values, not the values themselves
		if(fillRand == 3) //Probability 0.25 of adding a 4, and 0.75 of adding a 2
			board = fillSquare(board, fillIndex, toFill[1]);
		else
			board = fillSquare(board, fillIndex, toFill[0]);
		
		//printf("Updated board:\n");
		//printf("%ld\n",board);
		//printBoard_long(board);
		return board;
	}
	else return board;
}

void printBoard_long(unsigned long board)
{
	int newBoard[_size][_size];
	newPowerBoard(board, newBoard);
	printBoard_array(newBoard);
}

void printBoard_array(int board[][_size])
{
	for(int i = 0; i < _size; i++)
	{
		for(int j = 0; j < _size; j++)
		{
			printf("%4d",board[i][j]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n");
}

static inline unsigned long transpose(unsigned long board)
{
    return transpose_table[0][board & 0xffff] + transpose_table[1][(board >> 16) & 0xffff] + transpose_table[2][(board >> 32) & 0xffff] + transpose_table[3][(board >> 48) & 0xffff];
}

/*void expMoveCompress(unsigned long board, float expScore[], float currProb, int currDepth)
{
    expanded += 4;
    long newBoards[4];
    int heuristics[4];
    int empty;
    expanded += 4;
    for(int i = 0; i < _possMoves; i++)
    {
        newBoards[i] = doMove(board, allMoves[i]);
        
        empty = numEmpty(newBoards[i]);
        for(int fillIndex = 0; fillIndex < 16; fillIndex++)
        {
            
        }
    }
}*/

void expMove(unsigned long board, float expScore[], float currProb, int currDepth)
{
	for(int i = 0; i < _possMoves; i++)
	{
        expanded++;
		unsigned long newBoard = doMove(board, allMoves[i]);
        if(newBoard != board)
        {
            int empty = numEmpty(newBoard);
            if(empty > 0)
            {
                float newProb = currProb/empty;
                if(newProb < _minProb || currDepth == _maxDepth)
                    expScore[i] = (scoreBoard(newBoard) - currScore + heuristic(newBoard));
                else
                {
                    for(int fillIndex = 0; fillIndex < 16; fillIndex++)
                    {
                        if(((newBoard >> (fillIndex << 2)) & 0x0f) == 0)
                        {
                            for(int fillOption = 0; fillOption < fillOptions; fillOption++)
                            {
                                unsigned long boardFill = newBoard + ((long)toFill[fillOption] << (fillIndex << 2));
                                newProb = currProb*prob[fillOption]/empty;
                                expScore[i] += expMove_low(boardFill, newProb, currDepth + 1)*prob[fillOption]/empty;
                            }
                        }
                    }
                }
            }
        }
	}
}

//We don't have to keep track of which move is the best for further iterations, but only the highest value.
float expMove_low(unsigned long board, float currProb, int currDepth)
{
	float newMoveVal, maxMoveVal = -1;
	for(int i = 0; i < _possMoves; i++)
	{
        expanded++;
		unsigned long newBoard = doMove(board, allMoves[i]);
        if(newBoard != board)
        {
            int empty = numEmpty(newBoard);
            if(empty > 0)
            {
                if(currProb < (_minProb*empty) || currDepth == _maxDepth)
                    newMoveVal = scoreBoard(newBoard) - currScore + heuristic(newBoard);
                else
                {
                    newMoveVal = 0;
                    for(int fillIndex = 0; fillIndex < 16; fillIndex++)
                    {
                        if(((newBoard >> (fillIndex << 2)) & 0x0f) == 0)
                        {
                            for(int fillOption = 0; fillOption < fillOptions; fillOption++)
                            {
                                unsigned long boardFill = newBoard + ((long)toFill[fillOption] << (fillIndex << 2));
                                float newProb = currProb*prob[fillOption]/empty;
                                newMoveVal += expMove_low(boardFill, newProb, currDepth + 1)*prob[fillOption]/empty;
                            }
                        }
                    }
                }
                if(newMoveVal > maxMoveVal)
                    maxMoveVal = newMoveVal;
            }
        }
	}
    return maxMoveVal;
}

float heuristic(unsigned long board)
{
	float heurScore = 0.0;
    long flipBoard = transpose(board);

    heurScore += row_heur_table[board & 0xffff];
    heurScore += row_heur_table[flipBoard & 0xffff];
    
    heurScore += row_heur_table[(board >> 16) & 0xffff];
    heurScore += row_heur_table[(flipBoard >> 16) & 0xffff];
    
    heurScore += row_heur_table[(board >> 32) & 0xffff];
    heurScore += row_heur_table[(flipBoard >> 32) & 0xffff];
    
    heurScore += row_heur_table[(board >> 48) & 0xffff];
    heurScore += row_heur_table[(flipBoard >> 48) & 0xffff];
	
	//Add if the highest value is in a corner.
    /*float heurMaxCorner = 2000.0;
    int maxCol[4];
    maxCol[3] = max_table[board & 0xffff];
    maxCol[2] = max_table[(board >> 16) & 0xffff];
    maxCol[1] = max_table[(board >> 32) & 0xffff];
    maxCol[0] = max_table[(board >> 48) & 0xffff];

    int winners = (int)(((board >> (4*(3 - maxCol[0]))) & 0x0f)
                + (((board >> (16 + 4*(3 - maxCol[1]))) & 0x0f) << 4)
                + (((board >> (32 + 4*(3 - maxCol[2]))) & 0x0f) << 8)
                + (((board >> (48 + 4*(3 - maxCol[3]))) & 0x0f) << 12));
    int row = max_table[winners];
    int col = maxCol[row];
    //printf("(col1, col2, col3, col4): (%d, %d, %d, %d)\n", maxCol[0], maxCol[1], maxCol[2], maxCol[3]);
    //printf("Largest vals: ");
    //printf("(row, col): (%d, %d)\n",row,col);
    //printf("(%d, %d), (%d, %d)\n", row1, col1, row2, col2);
	if((row == 0 || row == _size - 1) && (col == 0 || col == _size - 1))
        heurScore += heurMaxCorner;*/
    
	return heurScore;
}

unsigned long doMove(unsigned long prevBoard, char move)
{
	unsigned long newBoard = 0;
	int lineVal;
	if(move == 'U')
	{
        unsigned long flipBoard = transpose(prevBoard);
        long newFlipBoard = 0;
		for(int i = 0; i < _size; i++)
		{
			lineVal = (flipBoard >> (i << 4)) & 0xffff;
			newFlipBoard |= ((unsigned long)row_lshift_table[lineVal]) << (i << 4);
		}
        newBoard = transpose(newFlipBoard);
	}
	else if(move == 'D')
	{
        unsigned long flipBoard = transpose(prevBoard);
        long newFlipBoard = 0;
		for(int i = 0; i < _size; i++)
		{
			lineVal = (flipBoard >> (i << 4)) & 0xffff;
			newFlipBoard |= ((unsigned long)row_rshift_table[lineVal]) << (i << 4);
		}
        newBoard = transpose(newFlipBoard);
	}
	else if(move == 'L')
	{
		for(int i = 0; i < _size; i++)
		{
			lineVal = (prevBoard >> (i << 4)) & 0xffff;
			newBoard |= ((unsigned long)row_lshift_table[lineVal]) << (i << 4);
		}
	}
	else if(move == 'R')
	{
		for(int i = 0; i < _size; i++)
		{
			lineVal = (prevBoard >> (i << 4)) & 0xffff;
			newBoard |= ((unsigned long)row_rshift_table[lineVal]) << (i << 4);
		}
	}
	return newBoard;
}

