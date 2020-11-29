#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@Time    : 2020-09-03 @IOP
@Author  : cndaqiang
@Blog    : cndaqiang.github.io
@File    : 读入含时xyz或POSCAR处理rdf
"""

import sys
import os
import numpy as np

#统计的最短最远距离
minr=0 #1e-15 #avoid rdf/0
maxr=10.0
dr=0.1
#因为rdf:g(r)=dN/(4*pi*dr*r^2*\rho_N), 一般仅在r=a_0时取最大值,此时dr就影响计算结果,dr翻倍就减半, dr越大越小,但是积分不变的
#
ngrid=int(np.ceil(maxr-minr)/dr)
#dr=(maxr-minr)/ngrid
ngrid=ngrid+2
grid=np.arange(ngrid)*dr+minr

#xyz文件时需要
#xyz格式类型,默认beforenum+1行是总原子数,header+1行开始是原子结构
#有的程序输出的格式不规范,在这里定义
beforenum=0 #/3 #beforenum+1行是总原子数
header=2    #/9 #header+1行开始是原子结构
#原胞
cell=np.zeros([3,3])
cell[0]=np.array([12,0,0]) #a
cell[1]=np.array([0,12,0]) #b
cell[2]=np.array([0,0,12]) #c

#-----Input File
if len(sys.argv) == 1:
    print("Usage: "+str(sys.argv[0])+" [ xxx.xyz | xxx.vasp ] [ none | a | a b c]")
    exit()
if len(sys.argv) > 1:
    inputfile = str(sys.argv[1])
else:
    inputfile = "./292.xyz"
if os.path.exists(inputfile):
    print("read from",inputfile)
else:
    print("EXIT: There are not "+inputfile)
    exit()
if len(sys.argv) == 3:  # 只输入一个晶格参数,视为cubic
    cell[:,:]=0
    for i in np.arange(3): cell[i,i]=float(sys.argv[2])
if len(sys.argv) == 5: #输入三个晶格参数,视为长方形 a b  c晶胞
    cell[:,:]=0
    for i in np.arange(3): cell[i,i]=float(sys.argv[2+i])

#-----


#----
def readxyz(xyzfile):
    f=open(xyzfile,'r')
    for i in np.arange(beforenum): f.readline()#原子数前有几行
    natoms=int(f.readline())#原子数
    ierror=f.seek(0,0)
    rownum=len(f.readlines())#总行数
    if rownum%(natoms+header) == 0:
        nstep=int(rownum/(natoms+header))  #每个xyz单元,2个
    else:
        print("rownum can't be split by natoms + header",rownum,natoms+header)
        nstep=int(rownum/(natoms+header))  #每个xyz单元,2个
    #nstep=2 也可以在这里限制总的nstep步数
    #收集元素种类和数量
    ntyp=np.array([]).astype(np.str)
    nat=np.array([])
    ierror=f.seek(0,0)
    for i in np.arange(header): f.readline()
    for i in np.arange(natoms):
        label=f.readline().split()[0]
        if np.where(ntyp==label)[0].size == 0:
           ntyp =   np.append(ntyp,label)
           nat  =   np.append(nat,1)
        else:
            nat[np.where(ntyp==label)[0][0]]=nat[np.where(ntyp==label)[0][0]]+1
    nat=nat.astype(np.int)
    #按照原子数排序,画图的顺序和中心原子
    ntyp=ntyp[np.argsort(nat)]
    nat=nat[np.argsort(nat)]

    natmax=int(np.max(nat))
    totaltype=ntyp.size
    #======
    xyz=np.zeros([nstep,totaltype,natmax,3]) #步数,元素,原子编号,xyz
    ierror=f.seek(0,0)
    for istep in np.arange(nstep):
        index=np.zeros(totaltype)
        for i in np.arange(header): f.readline()
        for iatom in np.arange(natoms):
            line=f.readline()
            label=line.split()[0]
            itype=np.where(ntyp==label)[0][0]
            for idir in np.arange(3):
                xyz[istep,itype,int(index[itype]),idir]=float(line.split()[idir+1])
            index[itype]=index[itype]+1
    f.close()
    return nstep,ntyp,nat,xyz
def readposcar(poscar):
    f=open(poscar,'r')
    f.readline()
    alat=np.float(f.readline().split()[0])
    cell=np.zeros([3,3])
    for i in np.arange(3):
        cell[i,0:3] = [ float(x) for x in f.readline().split()[0:3] ] 
    cell=cell*alat
    #
    nstep=1
    #元素
    ntyp=np.array(f.readline().split())
    nat=np.array(f.readline().split())
    totaltype=min(ntyp.size,nat.size)
    ntyp=ntyp[0:totaltype].copy()
    nat=nat[0:totaltype].copy().astype(np.int)
    #坐标类型
    xyztype=f.readline().split()[0]
    #读坐标
    natmax=int(np.max(nat))
    #坐标
    xyz=np.zeros([nstep,totaltype,natmax,3])
    for istep in np.arange(nstep):
        for itype in np.arange(totaltype):
            for i in np.arange(nat[itype]):
                xyz[istep,itype,i,0:3]=[ float(x) for x in f.readline().split()[0:3] ] 
                #
                #    xyz[istep,itype,i,0:3]=[ (xyz[istep,itype,i,0:3]*cell[:,d]).sum() for d in np.arange(3) ]


    #计算坐标
    if xyztype[0] == 'D' or xyztype[0] == 'd':
        Dxyz=xyz.copy()
        for d in np.arange(3):
            xyz[:,:,:,d]=(Dxyz*cell[:,d]).sum(axis=3)
        #for d in np.arange(3):
        #    xyz=xyz[:,:,:,0:3]*cell[0:3,d]
    else:
        xyz=xyz*alat
    return nstep,ntyp,nat,xyz,cell

def move2cell(nstep,ntyp,nat,xyz,cell):
    for istep in np.arange(nstep):
        for i in np.arange(ntyp.size):
            for j in np.arange(nat[i]):
                #下面这套算法不是严格的，暂时没有好的方法
                while xyz[istep,i,j,0] > cell[0,0]:
                    #print("i,j",i,j,"addcell 0")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]-cell[0,:]
                while xyz[istep,i,j,0] < 0:
                    #print("i,j",i,j,"delcell 0")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]+cell[0,:]
                while xyz[istep,i,j,1] > cell[1,1]:
                    #print("i,j",i,j,"addcell 1")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]-cell[1,:]
                while xyz[istep,i,j,1] < 0:
                    #print("i,j",i,j,"delcell 1")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]+cell[1,:]
                while xyz[istep,i,j,2] > cell[2,2]:
                    #print("i,j",i,j,"addcell 2")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]-cell[2,:]
                while xyz[istep,i,j,2] < 0:
                    #print("i,j",i,j,"delcell 2")
                    xyz[istep,i,j,:]=xyz[istep,i,j,:]+cell[2,:] 
def calV(cell):
    V=0
    V=V+(cell[0,1]*cell[1,2]-cell[0,2]*cell[1,1])*cell[2,0] 
    V=V+(cell[0,2]*cell[1,0]-cell[0,0]*cell[1,2])*cell[2,1]
    V=V+(cell[0,0]*cell[1,1]-cell[0,1]*cell[1,0])*cell[2,2]
    return V


def savedata(A,B,filename):
    f=open(filename,'w')
    for i in np.arange(A.size):
        f.write(str(A[i])+"\t"+str(B[i])+'\n')
    f.close()
#

def calrdf(A_in,B_in):
    A=A_in.copy()
    B=B_in.copy()
    rdf=np.zeros(ngrid)
    for i in np.arange(A.shape[0]):
        dis=np.sqrt(np.square(B-A[i]).sum(axis=1))
        choose=np.floor(  ( dis[(dis >= max(minr, 1e-5) )&(dis < maxr+10*dr )] - minr )/dr  ).astype(int)
        #忽略本身原子
        #向下取整,处以体积4/3pi(r^3 -(r-dr)^3)
        choose=choose[choose<ngrid]
        for uniq in np.unique(choose):
            rdf[uniq]=rdf[uniq]+choose[(choose==uniq)].size
    return rdf*1.0/A.shape[0] #算一个A原子就够了,算多个A原子进行平均
def rdf(nstep,ntyp,nat,xyz,cell):
    rdftype=int(np.arange(ntyp.size+1).sum())#几种rdf
    timerdf=np.zeros([nstep,rdftype,ngrid])#含时rdf
    rdflable=np.zeros(rdftype)#rdftype种键长
    rdflable=rdflable.astype(np.str)
    irdftype=0
    supercell=np.zeros(3)
    for i in np.arange(3):  supercell[i]=max(np.around(maxr/cell[i,i]),2) #超胞数量，至少3x3x3
    supercell.astype(np.int)
    print("Super cell",supercell)
    for istep in np.arange(nstep):
        for left in np.arange(ntyp.size):
            lindex=left#A元素编号
            A=xyz[istep,lindex,0:nat[lindex]].copy()#原胞中的A
            for right in np.arange(ntyp.size-left):
                rindex=lindex+right#B元素编号
                rdflable[irdftype]=ntyp[lindex]+"-"+ntyp[rindex]
                #下面计算AB两种原子的rdf
                for x in np.arange(-supercell[0],1+supercell[0]):
                    for y in np.arange(-supercell[1],1+supercell[1]):
                        for z in np.arange(-supercell[2],1+supercell[2]):
                            B=xyz[istep,rindex,0:nat[rindex],:]+x*cell[0]+y*cell[1]+z*cell[2]#超胞中的B
                            timerdf[istep,irdftype,0:ngrid] = timerdf[istep,irdftype,0:ngrid] + calrdf(A,B)
                rho=1.0*nat[rindex]/calV(cell)
                timerdf[istep,irdftype,0:ngrid]=timerdf[istep,irdftype,0:ngrid]/(rho*4.0/3.0*np.pi*dr*(3*grid[:]**2+3*grid[:]*dr+dr**2)) #dN/(4*pi*rho_N,r^2*dr)
                #这里使用体积4/3*pi[ (r+dr)^3-r^3 ]=4/3*pi[ dr^3 +3r^2dr+3r*dr^2]
                irdftype=irdftype+1
    print(rdflable)
    return rdftype,rdflable,timerdf #种类数,标签,数据timerdf[istep,irdftype,igrid]

if inputfile[-3:] == 'xyz':
    nstep,ntyp,nat,xyz=readxyz(inputfile)
    cell=cell.copy()#使用输入的cell
    print("XYZ: set cell to ")
    print('a:',cell[0])
    print('b:',cell[1])
    print('c:',cell[2])
else: #POSCAR
    del cell#删除为xyz设置的cell
    nstep,ntyp,nat,xyz,cell=readposcar(inputfile)

#move to unit cell
move2cell(nstep,ntyp,nat,xyz,cell)
print("nstep,ntyp,nat",nstep,ntyp,nat)

rdftype,rdflable,timerdf=rdf(nstep,ntyp,nat,xyz,cell)

avTrdf=timerdf.mean(axis=0)


#保存数据
for irdftype in np.arange(rdftype):
    rdffile=inputfile+"."+rdflable[irdftype]+'.averagerdf.dat'
    print('Save rdf to '+rdffile)
    savedata(np.around(grid,3),np.around(avTrdf[irdftype],6),rdffile)#分别保留3,5位有效数字

if False: #输出到屏幕
    for istep in np.arange(nstep):
        print("r,\t",rdflable,"\t,rdf")
        for igrid in np.arange(ngrid):#跳过初始的0
            print(np.around(grid[igrid],2),"\t",np.around(timerdf[istep,:,igrid],2))

#下面为画图脚本

if True: #画图
    from matplotlib import pyplot as plt
    #
    row=max(int(1.0*rdftype/np.sqrt(rdftype)),2) #向下取整
    col=max(int(np.ceil(rdftype/row)),2)
    fig, axs = plt.subplots(row,col,sharex=True,sharey=False,figsize=(8*col,8*row))
    i=0
    for x in np.arange(row):
        for y in np.arange(col):
            if i  < rdftype:#画单个键长
                axs[x,y].plot(grid,avTrdf[i],label=rdflable[i])
                axs[x,y].hlines(1.0, minr,maxr)
                axs[x,y].legend()
                i=i+1
            elif i  == rdftype :#如果还有空间就画所有键长
                for j in np.arange(rdftype):
                    axs[x,y].plot(grid,avTrdf[j],label=rdflable[j])
                    axs[x,y].hlines(1.0, minr,maxr)
                    axs[x,y].legend()
                    i=i+1
                break
            else: 
                break            
    plt.xlim(minr,maxr)
    plt.ylim(0,None)
    #plt.set_ylim(0,None)
    
    filename=inputfile+'.averagerdf.png'
    print('Save png to '+filename)
    plt.savefig(filename,dpi=100)