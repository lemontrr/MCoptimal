import affine_check
import math
import time
import datetime
import os
import sys
import multiprocessing
from optparse import OptionParser

def imp(n,La_x,Lb_z,Sbox,Abox,NOT_gates,kappa):
    imps = []
    space = La_x[:]
    for i in range(len(Abox)):
        g = Abox[i]; flg = False
        range_set = [b for b in range(1,1<<(n+i)) if bin(b).count('1')<kappa] 
        for a_i in range(len(range_set)): 
            a = range_set[a_i]
            for b_i in range(a_i+1,len(range_set)):
                b = range_set[b_i]
                if space[a]&space[b] == g:
                    imps.append((a,b,0))
                    flg = True
                    break
            if flg: break
        if not flg:
            for a in range(1,1<<(n+i)):
                for b in range(a+1,1<<(n+i)):
                    if space[a]&space[b] == g:
                        imps.append((a,b,0))
                        flg = True
                        break
                if flg: break
        if not flg:
            for a in range(1,1<<(n+i)):
                for b in range(a+1,1<<(n+i)):
                    for c in range(1<<(n+i)):
                        if (space[a]&space[b])^space[c] == g:
                            imps.append((a,b,c))
                            flg = True
                            break
                    if flg: break
                if flg: break
        space += [f^g for f in space]
    S_imps = []
    for i in range(n):
        for a in range(len(space)):
            if space[a] == Lb_z[1<<i]:
                S_imps.append(a)
                break
        
    text = f'S = []\nfor X in range({1<<n}):\n\tx=[(X>>i)&1 for i in range({n})]\n\tg = [0]*{len(Abox)}; y = [0]*{n}\n\n'
    for i in range(len(Abox)):
        a,b,c = imps[i]
        
        temp = f'\tg[{i}] = ('
        bin_a = bin(a)[2:].zfill(n+i)[::-1]
        for j in range(n+i):
            if bin_a[j] == '1':
                if j < n:
                    temp += f'x[{j}]^'
                else:
                    temp += f'g[{j-n}]^'
        temp = temp[:-1] + ') & ('
        bin_b = bin(b)[2:].zfill(n+i)[::-1]
        for j in range(n+i):
            if bin_b[j] == '1':
                if j < n:
                    temp += f'x[{j}]^'
                else:
                    temp += f'g[{j-n}]^'
        temp = temp[:-1] + ') '
        bin_c = bin(c)[2:].zfill(n+i)[::-1]
        for j in range(n+i):
            if bin_c[j] == '1':
                if j < n:
                    temp += f'^ x[{j}] '
                else:
                    temp += f'^ g[{j-n}] '
        temp += '\n'
        text += temp
    if len(Abox) != 0:
        text += '\n'
    for i in range(n):
        a = S_imps[i]
        temp = f'\ty[{i}] ='
        bin_a = bin(a)[2:].zfill(n+len(Abox))[::-1]
        for j in range(n+len(Abox)):
            if bin_a[j] == '1':
                if j < n:
                    temp += f' x[{j}] ^'
                else:
                    temp += f' g[{j-n}] ^'
        if NOT_gates[i] == 1:
            temp += ' 1 ^'
        temp = temp[:-1] + '\n'
        text += temp
    temp = '\tS.append('
    for i in range(n):
        temp += f'(y[{i}]<<{i})|'
    temp = temp[:-1] + ')\n'
    text += temp
    text += f'\nprint(S)\nprint({Sbox})\nif tuple(S)==tuple({Sbox}):\n\tprint(True)\nelse:\n\tprint(False)\n'
    return text

def Calculate_Remained_Targets(current_node,target_candidates):
    targets = []
    for i in target_candidates:
        if i not in current_node:
            targets.append(i)
    return targets
    

def plus_3AND(Abox_pool,targets_in,kappa):
    st = time.time()
    targets = Calculate_Remained_Targets(Abox_pool,targets_in)
    size = len(Abox_pool)
    target_pool = set()
    for f in targets:
        target_pool.update([i^f for i in Abox_pool])
    
    already_pool = set()
    subreturn_pool = []
    
    if kappa == 999:
        for b1 in range(1,size>>1):
            for b2 in range(b1+1,size): 
                f = Abox_pool[b1]&Abox_pool[b2]
                space = Abox_pool + [g^f for g in Abox_pool]
                if len(set(space)) < len(Abox_pool)*2: continue
                space.sort()
                tups = []; subspace = set(Abox_pool)
                for i in space:
                    if i not in subspace:
                        tups.append(i)
                        subspace.update([j^i for j in subspace])
                tups = tuple(tups)
                if tups not in already_pool:
                    already_pool.add(tups)
                    subreturn_pool.append(f)
    else:    
        range_set = [b for b in range(1,size) if bin(b).count('1')<kappa] 
        for b1_i in range(len(range_set)): 
            b1 = range_set[b1_i]
            for b2_i in range(b1_i+1,len(range_set)):
                b2 = range_set[b2_i] 
                f = Abox_pool[b1]&Abox_pool[b2]
                space = Abox_pool + [g^f for g in Abox_pool]
                if len(set(space)) < len(Abox_pool)*2: continue
                space.sort()
                tups = []; subspace = set(Abox_pool)
                for i in space:
                    if i not in subspace:
                        tups.append(i)
                        subspace.update([j^i for j in subspace])
                tups = tuple(tups)
                if tups not in already_pool:
                    already_pool.add(tups)
                    subreturn_pool.append(f)
    already_pool = set()
    return_pool = []

    subsubreturn_pool = []
    nmasd = 0
    if kappa == 999:
        for f in subreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [_^f for _ in Abox_pool]
            new_target_pool = list(target_pool) + [_^f for _ in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1 in range(1,size):
                for b2 in range(b1+1,size<<1): 
                    g = fspace[b1]&fspace[b2]
                    fgspace = fspace + [g^_ for _ in fspace]
                    if len(set(fgspace)) < len(fspace)*2: continue
                    fgspace.sort()
                    tups = []; subspace = set(Abox_pool)
                    for i in fgspace:
                        if i not in subspace:
                            tups.append(i)
                            subspace.update([j^i for j in subspace])
                    tups = tuple(tups)
                    if tups not in already_pool:
                        already_pool.add(tups)
                        subsubreturn_pool.append((f,g))
        
    else:
        range_set = [b for b in range(1,size<<1) if bin(b).count('1')<kappa] 
        for f in subreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [_^f for _ in Abox_pool]
            new_target_pool = list(target_pool) + [_^f for _ in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1_i in range(len(range_set)): 
                b1 = range_set[b1_i]
                for b2_i in range(b1_i+1,len(range_set)):
                    b2 = range_set[b2_i] 
                    g = fspace[b1]&fspace[b2]
                    fgspace = fspace + [g^_ for _ in fspace]
                    if len(set(fgspace)) < len(fspace)*2: continue
                    fgspace.sort()
                    tups = []; subspace = set(Abox_pool)
                    for i in fgspace:
                        if i not in subspace:
                            tups.append(i)
                            subspace.update([j^i for j in subspace])
                    tups = tuple(tups)
                    if tups not in already_pool:
                        already_pool.add(tups)
                        subsubreturn_pool.append((f,g))
                        
    already_pool = set()
    return_pool = []
    
    if kappa == 999:
        nmasd = 0
        for f,g in subsubreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [f^_ for _ in Abox_pool]
            fgspace = fspace + [g^_ for _ in fspace]
            new_target_pool = list(target_pool) + [f^_ for _ in target_pool] + [g^_ for _ in target_pool] + [f^g^_ for _ in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1 in range(1,size<<1):
                for b2 in range(b1+1,size<<2): 
                    h = fgspace[b1]&fgspace[b2]
                    
                    if h in new_target_pool:
                        space = fgspace + [g^_ for _ in fgspace]
                        space.sort()
                        tups = []; subspace = set(Abox_pool)
                        for i in space:
                            if i not in subspace:
                                tups.append(i)
                                subspace.update([j^i for j in subspace])
                        tups = tuple(tups)
                        if tups not in already_pool:
                            already_pool.add(tups)
                            return_pool.append((f,g,h))
    else:
        nmasd = 0
        range_set = [b for b in range(1,size<<2) if bin(b).count('1')<kappa]
        for f,g in subsubreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [f^_ for _ in Abox_pool]
            fgspace = fspace + [g^_ for _ in fspace]
            new_target_pool = list(target_pool) + [f^_ for _ in target_pool] + [g^_ for _ in target_pool] + [f^g^_ for _ in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1_i in range(len(range_set)): 
                b1 = range_set[b1_i]
                for b2_i in range(b1_i+1,len(range_set)):
                    b2 = range_set[b2_i] 
                    h = fgspace[b1]&fgspace[b2]
                    
                    if h in new_target_pool:
                        space = fgspace + [g^_ for _ in fgspace]
                        space.sort()
                        tups = []; subspace = set(Abox_pool)
                        for i in space:
                            if i not in subspace:
                                tups.append(i)
                                subspace.update([j^i for j in subspace])
                        tups = tuple(tups)
                        if tups not in already_pool:
                            already_pool.add(tups)
                            return_pool.append((f,g,h))
                    
                            if len(target_pool)==int(len(targets_in)/2):
                                return return_pool
    return return_pool

def plus_2AND(Abox_pool,targets_in,kappa):
    st = time.time()
    targets = Calculate_Remained_Targets(Abox_pool,targets_in)
    size = len(Abox_pool)
    target_pool = set()
    for f in targets:
        target_pool.update([i^f for i in Abox_pool])
    
    already_pool = set()
    subreturn_pool = []
    
    if kappa == 999:
        for b1 in range(1,size>>1):
            for b2 in range(b1+1,size): 
                f = Abox_pool[b1]&Abox_pool[b2]
                space = Abox_pool + [g^f for g in Abox_pool]
                if len(set(space)) < len(Abox_pool)*2: continue
                space.sort()
                tups = []; subspace = set(Abox_pool)
                for i in space:
                    if i not in subspace:
                        tups.append(i)
                        subspace.update([j^i for j in subspace])
                tups = tuple(tups)
                if tups not in already_pool:
                    already_pool.update([tups])
                    subreturn_pool.append(f)
    else:
        range_set = [b for b in range(1,size) if bin(b).count('1')<kappa]
        for b1_i in range(len(range_set)): 
            b1 = range_set[b1_i]
            for b2_i in range(b1_i+1,len(range_set)):
                b2 = range_set[b2_i] 
                f = Abox_pool[b1]&Abox_pool[b2]
                space = Abox_pool + [g^f for g in Abox_pool]
                if len(set(space)) < len(Abox_pool)*2: continue
                space.sort()
                tups = []; subspace = set(Abox_pool)
                for i in space:
                    if i not in subspace:
                        tups.append(i)
                        subspace.update([j^i for j in subspace])
                tups = tuple(tups)
                if tups not in already_pool:
                    already_pool.update([tups])
                    subreturn_pool.append(f)
    
    already_pool = set()
    return_pool = []
    
    if kappa == 999:
        range_set = [b for b in range(1,size<<1) if bin(b).count('1')<kappa] 
        nmasd = 0
        for f in subreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [g^f for g in Abox_pool]
            new_target_pool = list(target_pool) + [g^f for g in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1 in range(1,size):
                for b2 in range(b1+1,size<<1):
                    g = fspace[b1]&fspace[b2]
                    
                    if g in new_target_pool:
                        space = fspace + [g^i for i in fspace]
                        space.sort()
                        tups = []; subspace = set(Abox_pool)
                        for i in space:
                            if i not in subspace:
                                tups.append(i)
                                subspace.update([j^i for j in subspace])
                        tups = tuple(tups)
                        if tups not in already_pool:
                            already_pool.update([tups])
                            return_pool.append((f,g))
                    
                            if len(target_pool)==int(len(targets_in)/2):
                                return return_pool
    else:
        range_set = [b for b in range(1,size<<1) if bin(b).count('1')<kappa] 
        nmasd = 0
        for f in subreturn_pool:
            nmasd += 1
            fspace = Abox_pool + [g^f for g in Abox_pool]
            new_target_pool = list(target_pool) + [g^f for g in target_pool]
            new_target_pool = set(new_target_pool)
            
            for b1_i in range(len(range_set)): 
                b1 = range_set[b1_i]
                for b2_i in range(b1_i+1,len(range_set)):
                    b2 = range_set[b2_i] 
                    g = fspace[b1]&fspace[b2]
                    
                    if g in new_target_pool:
                        space = fspace + [g^i for i in fspace]
                        space.sort()
                        tups = []; subspace = set(Abox_pool)
                        for i in space:
                            if i not in subspace:
                                tups.append(i)
                                subspace.update([j^i for j in subspace])
                        tups = tuple(tups)
                        if tups not in already_pool:
                            already_pool.update([tups])
                            return_pool.append((f,g))
                    
                            if len(target_pool)==int(len(targets_in)/2):
                                return return_pool
    return return_pool

def plus_1AND(Abox_pool,targets_in,kappa):
    targets = Calculate_Remained_Targets(Abox_pool,targets_in)
    size = len(Abox_pool)
    target_pool = set()
    for f in targets:
        target_pool.update([i^f for i in Abox_pool])
    possible_pool = []
    
    if kappa == 999:
        for b1 in range(1,size>>1):
            for b2 in range(b1+1,size): 
                f = Abox_pool[b1]&Abox_pool[b2]
                if f in target_pool:
                    possible_pool.append(f)
    else:
        range_set = [b for b in range(1,size) if bin(b).count('1')<kappa] 
        for b1_i in range(len(range_set)): 
            b1 = range_set[b1_i]
            for b2_i in range(b1_i+1,len(range_set)):
                b2 = range_set[b2_i]
                f = Abox_pool[b1]&Abox_pool[b2]
                if f in target_pool:
                    possible_pool.append(f)

    already_pool = set()
    return_pool = []
    for f in possible_pool:
        space = Abox_pool + [g^f for g in Abox_pool]
        space.sort()
        tups = []; subspace = set(Abox_pool)
        for i in space:
            if i not in subspace:
                tups.append(i)
                subspace.update([j^i for j in subspace])
        tups = tuple(tups)
        if tups not in already_pool:
            already_pool.update([tups])
            return_pool.append(f)
    return return_pool

def A_star(n,Sbox,Sname,kappa=999, search_1AND_kappa = False, search_2AND = True, search_3AND = False, first_affine_search = True, excute_indexing = 'test',mpc = 8,folder='test'):
    start_time = datetime.datetime.now()
    st = time.time()
    NOT_gates = [(Sbox[0]>>i)&1 for i in range(n)]
    S = [Sbox[0]^Sbox[x] for x in range(1<<n)]
    size = 1<<n
    La_x = [0 for La in range(size)]
    Lb_z = [0 for Lb in range(size)]
    for La in range(size):
        for x in range(size):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
            Lb_z[La] |= (bin(La&S[x]).count('1')%2)<<x
    UB = [0,0,6,3+4,5+5,17+6,31+7,50+8,77+9,122+10,190+11]
    openList = [[[] for depth in range(n+2)] for f in range(UB[n])]
    openNum = [[0 for depth in range(n+2)] for f in range(UB[n])]
    openNu = [0 for f in range(UB[n])]
    openN = 0

    non_linears = [f for f in Lb_z if f not in La_x]
    
    L = int(math.log2(size-len(non_linears)))
    if L == n:
        text = f'# {tuple()}\n' + imp(n,La_x,Lb_z,Sbox,tuple(),NOT_gates,kappa)
        with open(f'./{folder}/{Sname}_{n}bit_{len(tuple())}ANDs_{time.time()-st:.2f}s.py','w') as fle:
            fle.write(text)
        return True


    if first_affine_search:
        D_info = dict(); all_D = dict()
        for mc_dim in [[1,2],[2,3],[2,4],[3,4],[3,5],[3,6],[4,5],[4,6],[4,7],[4,8]]:
            D = affine_check.Read_Data_set_1(n,[mc_dim])
            for i in D:
                D_info[i] = mc_dim
                all_D[i] = affine_check.Read_Data_Abox_set(n,i)
        
        if n == 3: min_MC,fAboxes = affine_check.First_find_affine_recipe_3bit(n,all_D,D_info,non_linears)
        elif n == 4: min_MC,fAboxes = affine_check.First_find_affine_recipe_4bit(n,all_D,D_info,non_linears)
        elif n == 5: min_MC,fAboxes = affine_check.First_find_affine_recipe_5bit(n,all_D,D_info,non_linears)
        elif n == 6: min_MC,fAboxes = affine_check.First_find_affine_recipe_6bit(n,all_D,D_info,non_linears)
        elif n == 7: min_MC,fAboxes = affine_check.First_find_affine_recipe_7bit(n,all_D,D_info,non_linears)
        elif n == 8: min_MC,fAboxes = affine_check.First_find_affine_recipe_8bit(n,all_D,D_info,non_linears)

        for func in fAboxes:
            for Abox in fAboxes[func]:
                f = len(Abox)+n-L-1
                d = 1
                
                openList[f][d].append(Abox)
                openNum[f][d] += 1
                openNu[f] += 1
                openN += 1
        f = min_MC+n-L-1; d = 1; bef_d = 0; openN_bef = openN
    else:
        f = n-L-1
        d = 0
        openList[f][d].append([])
        openNum[f][d] += 1
        openNu[f] += 1
        openN += 1
        min_MC = 1; bef_d = 0; openN_bef = openN
    pool = multiprocessing.Pool(processes=mpc)
    while openN > 0:
        if openNu[f] == 0:
            f += 1; d = n
        elif openNum[f][d] == 0: 
            d-=1
            if d == 0:
                d = n
        else:
            if bef_d < d:
                bef_d = d
            currentNode = openList[f][d][0]
            openList[f][d].remove(currentNode)
            openN -= 1
            openNu[f] -= 1
            openNum[f][d] -= 1
            if search_3AND and (len(currentNode)>1) and (currentNode[-1] == 0) and (currentNode[-2] == 0):
                ### MULTI CORE ###
                multi_nodes = [currentNode[:-2]]
                for i in range(len(openList[f][d])):
                    another_node = openList[f][d][i]
                    if (len(another_node)>1) and (another_node[-1] == 0) and (another_node[-2] == 0):
                        multi_nodes.append(another_node[:-2])
                        if len(multi_nodes) == mpc:
                            break
                for another_node in multi_nodes[1:]:
                    openList[f][d].remove(another_node+[0,0])
                    openN -= 1
                    openNu[f] -= 1
                    openNum[f][d] -= 1
                pc = len(multi_nodes)
                
                currentspace = [La_x[:] for i in range(pc)]
                for pc_i in range(pc):
                    for g in multi_nodes[pc_i]:
                        currentspace[pc_i] += [i^g for i in currentspace[pc_i]]
                
                results = [pool.apply_async(plus_3AND, (currentspace[pc_i], Lb_z, kappa)) for pc_i in range(pc)]
                output = [p.get() for p in results]
                
                New_fg = [[] for pc_i in range(pc)]
                for pc_i in range(pc):
                    New_fg[pc_i] += output[pc_i]
                                
                NewNodes = []
                for pc_i in range(pc):
                    NewNodes += [multi_nodes[pc_i] + [f,g] for f,g in New_fg[pc_i]]
                len_NewNodes = len(NewNodes)
                openList[f][d] = NewNodes + openList[f][d] 
                openN += len_NewNodes
                openNu[f] += len_NewNodes
                openNum[f][d] += len_NewNodes
                
            elif search_2AND and (len(currentNode)>0) and (currentNode[-1] == 0):
                ### MULTI CORE ###
                multi_nodes = [currentNode[:-1]]
                for i in range(len(openList[f][d])):
                    another_node = openList[f][d][i]
                    if len(another_node)+n-L-d != f:
                        multi_nodes.append(another_node[:-1])
                        if len(multi_nodes) == mpc:
                            break
                for another_node in multi_nodes[1:]:
                    openList[f][d].remove(another_node+[0])
                    openN -= 1
                    openNu[f] -= 1
                    openNum[f][d] -= 1
                pc = len(multi_nodes)
                
                currentspace = [La_x[:] for i in range(pc)]
                for pc_i in range(pc):
                    for g in multi_nodes[pc_i]:
                        currentspace[pc_i] += [i^g for i in currentspace[pc_i]]
                
                results = [pool.apply_async(plus_2AND, (currentspace[pc_i], Lb_z, kappa)) for pc_i in range(pc)]
                output = [p.get() for p in results]
                
                New_fg = [[] for pc_i in range(pc)]
                for pc_i in range(pc):
                    New_fg[pc_i] += output[pc_i]
                
                NewNodes = []
                for pc_i in range(pc):
                    NewNodes += [multi_nodes[pc_i] + [f,g] for f,g in New_fg[pc_i]]
                len_NewNodes = len(NewNodes)
                openList[f][d] = NewNodes + openList[f][d] 
                openN += len_NewNodes
                openNu[f] += len_NewNodes
                openNum[f][d] += len_NewNodes
                
                if search_3AND:                
                    for node in multi_nodes:
                        openList[f+1][d+1].append(node+[0,0])
                        openN += 1
                        openNu[f+1] += 1
                        openNum[f+1][d+1] += 1
                d = d+1
            elif d==n-L:
                text = f'# {currentNode}\n' + imp(n,La_x,Lb_z,Sbox,currentNode,NOT_gates,kappa)
                
                with open(f'./{folder}/{Sname}_{n}bit_{len(currentNode)}ANDs_{time.time()-st:.2f}s.py','w') as fle:
                    fle.write(text)
                print('finished_SUCCESS -',Sname,f'{time.time()-st:.2f}')
                with open(f'./{folder}/log_{excute_indexing}.txt','a') as fle:
                    fle.write(f'{Sname} ({len(currentNode)}ANDs) ::: {str(start_time)[:-4]} -> {str(datetime.datetime.now())[:-4]} (+{time.time()-st:.2f}s)\n')
                return True
            else:
                ### MULTI CORE ###
                multi_nodes = [currentNode]
                for i in range(len(openList[f][d])):
                    another_node = openList[f][d][i]
                    if another_node[-1] != 0:
                        multi_nodes.append(another_node)
                        if len(multi_nodes) == mpc:
                            break
                for another_node in multi_nodes[1:]:
                    openList[f][d].remove(another_node)
                    openN -= 1
                    openNu[f] -= 1
                    openNum[f][d] -= 1
                pc = len(multi_nodes)
                
                currentspace = [La_x[:] for i in range(pc)]
                for pc_i in range(pc):
                    for g in multi_nodes[pc_i]:
                        currentspace[pc_i] += [i^g for i in currentspace[pc_i]]
                if search_1AND_kappa:
                    results = [pool.apply_async(plus_1AND, (currentspace[pc_i], Lb_z, kappa)) for pc_i in range(pc)]
                else:
                    results = [pool.apply_async(plus_1AND, (currentspace[pc_i], Lb_z, 999)) for pc_i in range(pc)]
                output = [p.get() for p in results]
                
                New_g = [[] for pc_i in range(pc)]
                for pc_i in range(pc):
                    New_g[pc_i] += output[pc_i]
                    
                NewNodes = []
                for pc_i in range(pc):
                    NewNodes += [multi_nodes[pc_i] + [g] for g in New_g[pc_i]]
                len_NewNodes = len(NewNodes)
                openList[f][d+1] += NewNodes
                openN += len_NewNodes
                openNu[f] += len_NewNodes
                openNum[f][d+1] += len_NewNodes
                
                if search_2AND:
                    for node in multi_nodes:
                        openList[f+1][d+1].append(node+[0])
                        openN += 1
                        openNu[f+1] += 1
                        openNum[f+1][d+1] += 1
                d = d+1
    
    with open(f'./{folder}/log_{excute_indexing}.txt','a') as fle:
        fle.write(f'{Sname} - FAIL ::: {str(start_time)[:-4]} -> {str(datetime.datetime.now())[:-4]} (+{time.time()-st:.2f}s)\n')
    print('finished_FAIL -',Sname,f'{time.time()-st:.2f}')
    return False

def Sbox_formatting(S_text):
    values = S_text.replace('\n','').split(',')
    format_is = 'INT'
    for i in values:
        if 'x' in i:
            format_is = '0xHEX'
            continue
        for j in i:
            if j in ['a','b','c','d','e','f','A','B','C','D','E','F']:
                format_is = 'HEX'
    
    try:
        if format_is == 'INT':
            Sbox = list(map(lambda x:int(x),values))
        elif format_is == 'HEX':
            Sbox = list(map(lambda x:int('0x'+x,16),values))
        elif format_is == '0xHEX':
            Sbox = list(map(lambda x:int(x,16),values))
    except:
        Sbox = format_is

    return Sbox

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-S",dest="Sbox",help="S-box to implement (default identity)",action="store",default=f'{[i for i in range(32)]}'[1:-1])
    parser.add_option("-N",dest="Sname",help="S-box name (default test)",action="store",default='test')
    parser.add_option("-n",dest="size",help="Size of the S-box (default 5)",action="store",default=5)
    parser.add_option("-a",dest="op1",help="Use affine equivalence option",action="store_true",default=False)
    parser.add_option("-b",dest="k",help="Use kappa option",action="store",default=999)
    parser.add_option("-c",dest="kw",help="Use kappa option (without 1AND search)",action="store",default=999)
    parser.add_option("-M",dest="Multi",help="multi_threading",action="store",default=1)
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        print(parser.print_help())
        exit()
        
    if (int(options.k)<999) and (int(options.kw)<999):
        print("ERROR ::: Options -b and -c are incompatible.")
        exit()
        
    n = int(options.size)
    execute_list = []
    if '.txt' in options.Sbox:
        try:
            with open(options.Sbox,'r') as f:
                lines = f.readlines()
        except:
            print(f'ERROR ::: There is no file {options.Sbox}')
            exit()
        for i in range(len(lines)):
            line = lines[i]
            
            Sname,Sbox_text = line.split('=')
            Sbox = Sbox_formatting(Sbox_text)
            if type(Sbox) == str:
                print(f'ERROR ::: The file is in the wrong format at line {i+1}')
                print(f'ERROR ::: Its format is assumed to be {Sbox}')
                exit()
            execute_list.append([Sname,Sbox])
                
    else:
        Sbox = Sbox_formatting(options.Sbox)
        execute_list.append([options.Sname,Sbox])
    
    for Sname,Sbox in execute_list:
        if len(Sbox) != 1<<n:
            print('ERROR ::: The size and S-box are not compatible.')
            print(f'ERROR ::: The length of S-box is {len(Sbox)}, but you input as size n={n}')
            print(f'ERROR ::: name of S-box : {Sname}')
            print(f'ERROR ::: S-box : {Sbox}')
            exit()
        else:
            for i in range(len(Sbox)):
                if Sbox[i]>(1<<n):
                    print('ERROR ::: The size and S-box are not compatible.')
                    print(f'ERROR ::: There is an invalid value S[{i}]={Sbox[i]}')
                    print(f'ERROR ::: name of S-box : {Sname}')
                    print(f'ERROR ::: S-box : {Sbox}')
                    exit()
    
    filelist = os.listdir('./')
    file_index = 0
    while(True):
        if (file_index % 10) == 1:
            indexing = f'{file_index}st'
        elif (file_index % 10) == 2:
            indexing = f'{file_index}nd'
        elif (file_index % 10) == 3:
            indexing = f'{file_index}rd'
        else:
            indexing = f'{file_index}th'
        
        is_there = False
        for file_name in filelist:
            if f'result_{indexing}' in file_name:
                is_there = True
        if is_there:
            file_index += 1
        else:
            break
        
    if (file_index % 10) == 1:
        indexing = f'{file_index}st'
    elif (file_index % 10) == 2:
        indexing = f'{file_index}nd'
    elif (file_index % 10) == 3:
        indexing = f'{file_index}rd'
    else:
        indexing = f'{file_index}th'
    
    folder = f'result_{indexing}_{n}bit'
    if options.op1:
        folder += f'_a'
    if int(options.k) < 999:
        search_1AND_kappa = True
        kappa = int(options.k)
        folder += f'_b{int(options.k)}'
    elif int(options.kw) < 999:
        search_1AND_kappa = False
        kappa = int(options.kw)
        folder += f'_c{int(options.kw)}'
    else:
        kappa = 999
        search_1AND_kappa = False
        
    print(f'{indexing} start!')    
    os.mkdir(folder)
    for Sname,Sbox in execute_list:        
        A_star(n,Sbox,Sname,kappa=kappa,search_1AND_kappa=search_1AND_kappa,first_affine_search=options.op1,excute_indexing=indexing,mpc=int(options.Multi),folder=folder)
    