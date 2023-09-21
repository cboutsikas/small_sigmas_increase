from scipy import io
import numpy as np
from scipy import sparse
import os
from matplotlib import pyplot as plt

# this python script creates plots (and store them as .pdfs).
# you need to give manually, the name, which is the string
# returned in ComputeAndStore() function, the exponent (10 base)
# of the largest singular value s_max, and the 
# exponent (10 base) of the smallest singular value s_min.
# Finally, you need to spesify the string prec; If prec == "single",
# we store a .pdf figure for the computed singular values in single 
# separately (among with the fixed singular values and computed in double).
# If prec == "half", we store a .pdf figure for the computed singular values in half,
# (among with the fixed singular values and computed in double). Nevertheless,
# we always store a .pdf showing the computed singular values for each precision,
# and for the fixed ones. The corresponding figures are being stores into
# the /figures folder.

def print_sigmas(name,s_max,s_min,prec):
    #os.chdir("data")
    M = io.mmread("data\M_"+name+".mtx")
    M = sparse.coo_matrix.toarray(M)
    x = np.arange(256)+1
    os.chdir("figures")
    # print for all the precisions
    fig2, ax = plt.subplots(1,3)
    fig2.tight_layout()
    s = np.arange(s_min,s_max,step=2)
    plt.setp(ax, yticks = np.power(np.ones(np.size(s))*10,s))
    plt.yticks(np.power(np.ones(np.size(s))*10,s))
    plt.sca(ax[0]) 
    plt.yticks(np.power(np.ones(np.size(s))*10,s))
    ax[0].set_yscale('log')
    ax[0].grid(True)
    ax[0].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b', marker = 's',
                  label = 'Double')
    ax[0].scatter(x,M[:,0], facecolors = 'none', s = 60, edgecolors = 'm', marker = 'v',
                  label = 'Exact')
    ax[0].legend(loc="lower left")

    ax[1].grid(True)
    ax[1].set_yscale('log')
    ax[1].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b',marker="s",
                  label = 'Double')
    ax[1].scatter(x,M[:,2], facecolors = 'none', s = 60, edgecolors = 'tab:orange',
                  marker="*",label = 'Single')
    ax[1].legend(loc="lower left")

    ax[2].grid(True)
    ax[2].set_yscale('log')
    ax[2].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b',marker="s",
                  label = 'Double')
    ax[2].scatter(x,M[:,3], facecolors = 'none', s = 60, edgecolors = 'g',marker="o", 
                  label = 'Half')
    ax[2].legend(loc="lower left")
    #plt.savefig("figures\"+name+".pdf")
    plt.savefig(name+'_.pdf')

    if prec == "single":
        fig2, ax = plt.subplots(1,2)
        fig2.tight_layout()
        s = np.arange(s_min,s_max,step=2)
        plt.setp(ax, yticks = np.power(np.ones(np.size(s))*10,s))
        plt.yticks(np.power(np.ones(np.size(s))*10,s))
        plt.sca(ax[0]) 
        plt.yticks(np.power(np.ones(np.size(s))*10,s))
        ax[0].set_yscale('log')
        ax[0].grid(True)
        ax[0].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b', 
                      marker = 's',label = 'Double')
        ax[0].scatter(x,M[:,0], facecolors = 'none', s = 60, edgecolors = 'm',
                       marker = 'v',label = 'Exact')
        ax[0].legend(loc="lower left")

        ax[1].grid(True)
        ax[1].set_yscale('log')
        ax[1].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b',
                      marker="s",label = 'Double')
        ax[1].scatter(x,M[:,2], facecolors = 'none', s = 60, edgecolors = 'tab:orange',
                    marker="*",label = 'Single')
        ax[1].legend(loc="lower left")
        plt.savefig(name+'_single.pdf')

    if prec == "half":
        fig2, ax = plt.subplots(1,2)
        fig2.tight_layout()
        s = np.arange(-6,2,step=2)
        plt.setp(ax, yticks = np.power(np.ones(np.size(s))*10,s))
        plt.yticks(np.power(np.ones(np.size(s))*10,s))
        plt.sca(ax[0]) 
        plt.yticks(np.power(np.ones(np.size(s))*10,s))
        ax[0].set_yscale('log')
        ax[0].grid(True)
        ax[0].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b', 
                      marker = 's', label = 'Double')
        ax[0].scatter(x,M[:,0], facecolors = 'none', s = 60, edgecolors = 'm', 
                      marker = 'v', label = 'Exact')
        ax[0].legend(loc="lower left")

        ax[1].grid(True)
        ax[1].set_yscale('log')
        ax[1].scatter(x,M[:,1], facecolors = 'none', s = 60, edgecolors = 'b',
                      marker="s",label = 'Double')
        ax[1].scatter(x,M[:,3], facecolors = 'none', s = 60, edgecolors = 'g',
                      marker="*",label = 'Half')
        ax[1].legend(loc="lower left")
        plt.savefig(name+'_half.pdf')

    os.chdir("../")


# you need to manually tune this parameter
print_sigmas("2_228_28_2_2_2_2_-4",2,-4,"half")

