#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:54:14 2019

@author: pankhuri
"""

#TSSN Program
##Initial code is for STS and TST switches


##Output-  Error messages
#         complexity
#         complexity variation vs k,beta graph
#         complexity vairation vs beta for lee and jacobian


import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
import math

def s(switch, blocking, inlet, occp, beta,alpha):
        
    N = inlet;
    complexity = N*N;
    print(complexity)
    return complexity

def sss(switch, blocking, inlet, occp, beta,alpha):
    N =inlet;
    k = 7;
    n = math.sqrt(N/2)
    complexity = 2*N*k + k*math.pow(N/n,2);
    print("complexity",complexity)
    return complexity

def ts(switch, blocking, inlet, occp, beta,alpha):

    p = occp;
    N = inlet;
    Nx = N*N; #    Nx = Number of space stage cross points
    Nc = 24; #number of control words
    Nbx = N*Nc* 7; # Nbx = Number of memory bits for the space stage control store = N × (Number of control words) (number of bits per control word)
    Nbt = N*Nc*8 + N*Nc*5;1
    nonblk  = math.sqrt(N/2);
    k_nonblk = 2*nonblk - 1;
    beta_cal = k_nonblk/nonblk;
    complexity = Nx + (Nbx + Nbt)/100;
    K = 7;
    #Blocking probablity
    T = 128; #    T = number of time slots in a frame  
    
    if (blocking == 0):
        B = p/(N*math.pow(T,3)); 
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,k_nonblk)*math.pow(2-p,2*N-k_nonblk)/(math.factorial(k_nonblk)*math.factorial(2*N-k_nonblk));
    if (blocking == 1):
        k = 7;
        B = p/(N*math.pow(T,3)); 
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
        
    if (blocking == 2):
        k = 7;
        B = p/(N*math.pow(T,3)); 
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
    
    else:
        print("error");
    B = p/(N*math.pow(T,3));   
    #blocking probability of a three stage switch using Jacobeaus
    blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));

    print(complexity);
    plt.subplot(2,1,1)
    plt.plot(complexity, k_nonblk);

    return complexity

def sts(switch, blocking, inlet, occp, beta,alpha):

    p = occp;
    N = inlet;
    K = 7;    
    nonblk  = math.sqrt(N/2);
    k_nonblk = 2*nonblk - 1;
    beta_nonblk = k_nonblk/nonblk;
    beta_blk = K/N;
    C = 128  #message channels
    Nx = 2*K*N
    Nb = (2*K*C*math.log2(N)+ K*C*8 + K*C* math.log2(C))
    complexity = Nx + Nb/100;
    #k = 1; #(k = 1 in second stage)
    
    #Blocking probablity
    T = 128; #    T = number of time slots in a frame  
    if (blocking == 0):
        B = math.pow(float(1-math.pow((1-p/beta_nonblk),2)),k_nonblk);
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,k_nonblk)*math.pow(2-p,2*N-k_nonblk)/(math.factorial(k_nonblk)*math.factorial(2*N-k_nonblk));
    if (blocking == 1):
        B = math.pow(float(1-math.pow((1-p/beta_blk),2)),K);
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
        
    if (blocking == 2):
        B = math.pow(float(1-math.pow((1-p/beta_blk),2)),K);
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
    
    else:
        print("error");

    

    print("Complexity of implementation is", complexity);
    print("Blocking probability with lee's graph", B);
    print("Blocking probability with Jacobian method", blockingprob_jac );

    complexity_graph_plot = [];
    beta_blk_graph_plot = [];
    k_graph_plot = [];

    
    for k_graph in (1, 2, 3, 4, 5, 6, 7):
        beta_blk_graph = k_graph/N;
        Nx_graph = 2*k_graph*N;
        Nb_graph = (2*k_graph*C*math.log2(N)+ k_graph*C*8 + k_graph*C* math.log2(C));
        complexity_graph = Nx_graph + Nb_graph/100;
        B_graph = math.pow(float(1-math.pow((1-p/beta_blk_graph),2)),k_graph);
        blockingprob_jac_graph = math.pow(math.factorial(N),2)* math.pow(p,k_graph)*math.pow(2-p,2*N-k_graph)/(math.factorial(k_graph)*math.factorial(2*N-k_graph));
        complexity_graph_plot.append(complexity_graph);
        k_graph_plot.append(k_graph);
        beta_blk_graph_plot.append(beta_blk_graph);

     
    plt.subplot(2,1,1)
    
    plt.plot(complexity_graph_plot, k_graph_plot, '-b', label='k');
    plt.title('Complexity vs k grgraph')


    plt.subplots_adjust(hspace = 1)
    plt.subplot(2,1,2)    
    plt.plot(complexity_graph_plot, beta_blk_graph_plot, '-r', label='beta');
    plt.title('Complexity vs beta graph')
    
    return complexity

def tst(switch, blocking, inlet, occp, beta,alpha):

    p = occp;
    N = inlet;
    K = 7;    
    nonblk  = math.sqrt(N/2);
    k_nonblk = 2*nonblk - 1;
    beta_cal = k_nonblk/nonblk;
    betaa = K/N;
    C = 128  #message channels
    
    T = 128; #    T = number of time slots in a frame  
    
    #Blocking probablity
    if (blocking == 0):
        L = 2*T -1; #L = number of space slot of space switch.
        B = math.pow(float(1-math.pow((1-p*T/L),N)),L)
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
        Nx = N*N
        Nb = (N*L*math.log2(N)+ 2*N*T*8 + 2*N*L* math.log2(T))
    if (blocking == 1):
        L = 25;
        B = math.pow(float(1-math.pow((1-p*T/L),N)),L)
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
        Nx = N*N
        Nb = (N*L*math.log2(N)+ 2*N*T*8 + 2*N*L* math.log2(T))
    if (blocking == 2):
        L = 25;
        B = math.pow(float(1-math.pow((1-p*T/L),N)),L)
        blockingprob_jac = math.pow(math.factorial(N),2)* math.pow(p,K)*math.pow(2-p,2*N-K)/(math.factorial(K)*math.factorial(2*N-K));
        Nx = N*N
        Nb = (N*L*math.log2(N)+ 2*N*T*8 + 2*N*L* math.log2(T))
    else:
        print("error");    


    complexity = Nx + Nb/100;
    
    print("Complexity of implementation is", complexity,);
    print("Blocking probability with lee's graph", B);
    print("Blocking probability with Jacobian methods", blockingprob_jac );
    
    complexity_graph_plot = [];
    beta_blk_graph_plot = [];
    k_graph_plot = [];

    
    for k_graph in (1, 2, 3, 4, 5, 6, 7):
        beta_blk_graph = k_graph/N;
        Nx_graph = 2*k_graph*N;
        Nb_graph = (2*k_graph*C*math.log2(N)+ k_graph*C*8 + k_graph*C* math.log2(C));
        complexity_graph = Nx_graph + Nb_graph/100;
        B_graph = math.pow(float(1-math.pow((1-p/beta_blk_graph),2)),k_graph);
        blockingprob_jac_graph = math.pow(math.factorial(N),2)* math.pow(p,k_graph)*math.pow(2-p,2*N-k_graph)/(math.factorial(k_graph)*math.factorial(2*N-k_graph));
        complexity_graph_plot.append(complexity_graph);
        k_graph_plot.append(k_graph);
        beta_blk_graph_plot.append(beta_blk_graph);

     
    plt.subplot(2,1,1)
    
    plt.plot(complexity_graph_plot, k_graph_plot, '-b', label='k');
    plt.title('Complexity vs k grgraph')


    plt.subplots_adjust(hspace = 1)
    plt.subplot(2,1,2)    
    plt.plot(complexity_graph_plot, beta_blk_graph_plot, '-r', label='beta');
    plt.title('Complexity vs beta graph')

    return complexity

def main(): 
    print("Enter switch type: 0- S 1- TS; 2-STS 3- TST 4- SSS")
    switch = int(input())
    print("Enter blocking type 1- non-blocking; 2-blocking with lee's method 3- blocking with Jacobian method")
    blocking = int(input())
    print("number of inlets ")
    inlet = int(input())
    print("Enter occupancy probability ")
    occp = float(input())
    print("Enter beta: ")
    beta = float(input())
    print("Enter alpha")
    alpha = float(input())

    if (switch == 0):
        s(switch, blocking, inlet, occp, beta, alpha)     
    if (switch == 1):
        ts(switch, blocking, inlet, occp, beta, alpha)
    if (switch == 2):
        sts(switch, blocking, inlet, occp, beta, alpha)
    if (switch == 3):
        tst(switch, blocking, inlet, occp, beta, alpha)
    if (switch == 4):
        sss(switch, blocking, inlet, occp, beta, alpha)

   
# call main 
if __name__ == '__main__': 
	main() 
