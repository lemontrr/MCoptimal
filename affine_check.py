import re
import os
import math
def autocorrelation(n,B):
    size = 1<<n
    autoC_abs = []
    fx = B
    for s in range(size):
        fx_s = 0
        for x in range(size):
            fx_s |= ((fx>>(x^s))&1)<<x
        autoC_abs.append(abs(bin(fx_s^fx).count('1')-int(size/2)))
    return autoC_abs

def walshtransform(n,B):
    size = 1<<n
    La_x = [0 for La in range(size)]
    for La in range(size):
        for x in range(size):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    Walsh_abs = [abs(bin(La_x[La]^B).count('1')-int(size/2)) for La in range(size)]
    return Walsh_abs

def invariants(n,B):
    size = 1<<n
    Walsh_abs = walshtransform(n,B)
    autoC_abs = autocorrelation(n,B)
    Walsh_abs_dis = dict()
    autoC_abs_dis = dict()
    for x in Walsh_abs:
        if x not in Walsh_abs_dis:
            Walsh_abs_dis[x] = 1
        else:
            Walsh_abs_dis[x] += 1
    for x in autoC_abs:
        if x not in autoC_abs_dis:
            autoC_abs_dis[x] = 1
        else:
            autoC_abs_dis[x] += 1
    Walsh_abs_tup = []
    autoC_abs_tup = []
    for x in range(size):
        if x in Walsh_abs_dis:
            Walsh_abs_tup.append(x)
            Walsh_abs_tup.append(Walsh_abs_dis[x])
        if x in autoC_abs_dis:
            autoC_abs_tup.append(x)
            autoC_abs_tup.append(autoC_abs_dis[x])
    return tuple(Walsh_abs_tup),tuple(autoC_abs_tup)
def Derivation(n,a,B):
    size = 1<<n
    Daf = 0
    for x in range(size):
        Daf |= (((B>>x)&1)^((B>>(x^a))&1))<<x
    return Daf

def decompose(n,a,f):
    size = 1<<n
    f0 = 0; f1 = 0
    x_0 = 0; x_1 = 0
    for x in range(size):
        if bin(a&x).count('1')%2:
            f1 |= ((f>>x)&1)<<x_1
            x_1 += 1
        else:
            f0 |= ((f>>x)&1)<<x_0
            x_0 += 1
    return f0,f1

def First_find_affine_recipe_3bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
    
                        if min_MC > MC:
                            min_MC = MC
                    
                        Abox_set_g = affine_func[g]
                        for Abox_g in Abox_set_g:
                            Abox_f = []
                            for i in range(len(Abox_g)-1):
                                h = Abox_g[i]
                                table_h = [(h>>i)&1 for i in range(1<<n)]
                                hA = 0
                                for i in range(1<<n):
                                    xA = 0
                                    for j in range(n):
                                        xA ^= ((i>>j)&1)*A[j]
                                    hA |= table_h[xA]<<i
                                Abox_f.append(hA)
                            Abox_f.append(f)
                            if f not in targets_Aboxes:
                                targets_Aboxes[f] = [Abox_f]
                            else:
                                targets_Aboxes[f].append(Abox_f)
                        f_flag = True
                        break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def First_find_affine_recipe_4bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
                        for A[3] in A_candidate[3]:
                            t = table_g[A[3]]
                            for i in range(1,8):
                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^A[3]]<<i
                            g01234 = g0123 | (t<<8)
                            if g01234 not in fC_part[3]: continue
                            if min_MC > MC:
                                min_MC = MC
                        
                            Abox_set_g = affine_func[g]
                            for Abox_g in Abox_set_g:
                                Abox_f = []
                                for i in range(len(Abox_g)-1):
                                    h = Abox_g[i]
                                    # print(Abox_g)
                                    table_h = [(h>>i)&1 for i in range(1<<n)]
                                    hA = 0
                                    for i in range(1<<n):
                                        xA = 0
                                        for j in range(n):
                                            xA ^= ((i>>j)&1)*A[j]
                                        hA |= table_h[xA]<<i
                                    Abox_f.append(hA)
                                Abox_f.append(f)
                                if f not in targets_Aboxes:
                                    targets_Aboxes[f] = [Abox_f]
                                else:
                                    targets_Aboxes[f].append(Abox_f)
                            f_flag = True
                            break
                        if f_flag: break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def First_find_affine_recipe_5bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
                        for A[3] in A_candidate[3]:
                            t = table_g[A[3]]
                            for i in range(1,8):
                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^A[3]]<<i
                            g01234 = g0123 | (t<<8)
                            if g01234 not in fC_part[3]: continue
                        
                            for A[4] in A_candidate[4]:
                                t = table_g[A[4]]
                                for i in range(1,16):
                                    t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^A[4]]<<i
                                g012345 = g01234 | (t<<16)
                                if g012345 not in fC_part[4]: continue
                                if min_MC > MC:
                                    min_MC = MC
                        
                                Abox_set_g = affine_func[g]
                                for Abox_g in Abox_set_g:
                                    Abox_f = []
                                    for i in range(len(Abox_g)-1):
                                        h = Abox_g[i]
                                        table_h = [(h>>i)&1 for i in range(1<<n)]
                                        hA = 0
                                        for i in range(1<<n):
                                            xA = 0
                                            for j in range(n):
                                                xA ^= ((i>>j)&1)*A[j]
                                            hA |= table_h[xA]<<i
                                        Abox_f.append(hA)
                                    Abox_f.append(f)
                                    if f not in targets_Aboxes:
                                        targets_Aboxes[f] = [Abox_f]
                                    else:
                                        targets_Aboxes[f].append(Abox_f)
                                f_flag = True
                                break
                            if f_flag: break
                        if f_flag: break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def First_find_affine_recipe_6bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
                        
                        for A[3] in A_candidate[3]:
                            t = table_g[A[3]]
                            for i in range(1,8):
                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^A[3]]<<i
                            g01234 = g0123 | (t<<8)
                            if g01234 not in fC_part[3]: continue
                        
                            for A[4] in A_candidate[4]:
                                t = table_g[A[4]]
                                for i in range(1,16):
                                    t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^A[4]]<<i
                                g012345 = g01234 | (t<<16)
                                if g012345 not in fC_part[4]: continue
                        
                                for A[5] in A_candidate[5]:
                                    t = table_g[A[5]]
                                    for i in range(1,32):
                                        t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^A[5]]<<i
                                    g0123456 = g012345 | (t<<32)
                                    if g0123456 not in fC_part[5]: continue
                                    if min_MC > MC:
                                        min_MC = MC
                            
                                    Abox_set_g = affine_func[g]
                                    for Abox_g in Abox_set_g:
                                        Abox_f = []
                                        for i in range(len(Abox_g)-1):
                                            h = Abox_g[i]
                                            table_h = [(h>>i)&1 for i in range(1<<n)]
                                            hA = 0
                                            for i in range(1<<n):
                                                xA = 0
                                                for j in range(n):
                                                    xA ^= ((i>>j)&1)*A[j]
                                                hA |= table_h[xA]<<i
                                            Abox_f.append(hA)
                                        Abox_f.append(f)
                                        if f not in targets_Aboxes:
                                            targets_Aboxes[f] = [Abox_f]
                                        else:
                                            targets_Aboxes[f].append(Abox_f)
                                    f_flag = True
                                    break
                                if f_flag: break
                            if f_flag: break
                        if f_flag: break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def First_find_affine_recipe_7bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
                        
                        for A[3] in A_candidate[3]:
                            t = table_g[A[3]]
                            for i in range(1,8):
                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^A[3]]<<i
                            g01234 = g0123 | (t<<8)
                            if g01234 not in fC_part[3]: continue
                        
                            for A[4] in A_candidate[4]:
                                t = table_g[A[4]]
                                for i in range(1,16):
                                    t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^A[4]]<<i
                                g012345 = g01234 | (t<<16)
                                if g012345 not in fC_part[4]: continue
                        
                                for A[5] in A_candidate[5]:
                                    t = table_g[A[5]]
                                    for i in range(1,32):
                                        t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^A[5]]<<i
                                    g0123456 = g012345 | (t<<32)
                                    if g0123456 not in fC_part[5]: continue
                            
                            
                                    for A[6] in A_candidate[6]:
                                        t = table_g[A[6]]
                                        for i in range(1,64):
                                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^(((i>>5)&1)*A[5])^A[6]]<<i
                                        g01234567 = g0123456 | (t<<64)
                                        if g01234567 not in fC_part[6]: continue
                                        if min_MC > MC:
                                            min_MC = MC
                                
                                        Abox_set_g = affine_func[g]
                                        for Abox_g in Abox_set_g:
                                            Abox_f = []
                                            for i in range(len(Abox_g)-1):
                                                h = Abox_g[i]
                                                table_h = [(h>>i)&1 for i in range(1<<n)]
                                                hA = 0
                                                for i in range(1<<n):
                                                    xA = 0
                                                    for j in range(n):
                                                        xA ^= ((i>>j)&1)*A[j]
                                                    hA |= table_h[xA]<<i
                                                Abox_f.append(hA)
                                            Abox_f.append(f)
                                            if f not in targets_Aboxes:
                                                targets_Aboxes[f] = [Abox_f]
                                            else:
                                                targets_Aboxes[f].append(Abox_f)
                                        f_flag = True
                                        break
                                    if f_flag: break
                                if f_flag: break
                            if f_flag: break
                        if f_flag: break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def First_find_affine_recipe_8bit(n,affine_func,affine_func_info,target_func):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    targets_Aboxes = dict(); min_MC = 5
    for f in target_func:
        f_flag = False
        inv_f = invariants(n,f)
        for g in affine_func:
            MC,dim = affine_func_info[g]
            if min_MC <= MC: continue
            
            inv_g = invariants(n,g)
            if inv_f != inv_g: continue
            A_candidate = [[j for j in range(1<<dim)] for i in range(n)]
            
            autoC_f = autocorrelation(n,f)
            autoC_g = autocorrelation(n,g)
            for i in range(n):
                e_i = 1<<i
                A_i = []
                for v in A_candidate[i]:
                    if autoC_g[v] == autoC_f[e_i]:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                Df = Derivation(n,e_i,f)
                inv_Df = invariants(n,Df)
                A_i = []
                for v in A_candidate[i]:
                    Dg = Derivation(n,v,g)
                    inv_Dg = invariants(n,Dg)
                    if inv_Df == inv_Dg:
                        A_i.append(v)
                A_candidate[i] = A_i

            for i in range(n):
                e_i = 1<<i
                f_ei = f ^ (1<<e_i)
                inv_f_ei = invariants(n,f_ei)
                A_i = []
                for v in A_candidate[i]:
                    g_v = g ^ (1<<v)
                    inv_g_v = invariants(n,g_v)
                    if inv_f_ei == inv_g_v:
                        A_i.append(v)
                A_candidate[i] = A_i
                
                ########## A_candidate ##########
            
            fC = [i^f for i in La_x]
            fC_part = [[i&((1<<(1<<j))-1) for i in fC] for j in range(1,n+1)]
            table_g = [(g>>i)&1 for i in range(1<<n)]
            A = [0]*n; g0 = 0
            for A[0] in A_candidate[0]:
                t = table_g[A[0]]&1
                g01 = g0 | (t<<1)
                if g01 not in fC_part[0]: continue
                
                for A[1] in A_candidate[1]:
                    t = table_g[A[1]]
                    t |= table_g[A[0]^A[1]]<<1
                    g012 = g01 | (t<<2)
                    if g012 not in fC_part[1]: continue
                    
                    for A[2] in A_candidate[2]:
                        t = table_g[A[2]]
                        for i in range(1,4):
                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^A[2]]<<i
                        g0123 = g012 | (t<<4)
                        if g0123 not in fC_part[2]: continue
                        
                        for A[3] in A_candidate[3]:
                            t = table_g[A[3]]
                            for i in range(1,8):
                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^A[3]]<<i
                            g01234 = g0123 | (t<<8)
                            if g01234 not in fC_part[3]: continue
                        
                            for A[4] in A_candidate[4]:
                                t = table_g[A[4]]
                                for i in range(1,16):
                                    t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^A[4]]<<i
                                g012345 = g01234 | (t<<16)
                                if g012345 not in fC_part[4]: continue
                        
                                for A[5] in A_candidate[5]:
                                    t = table_g[A[5]]
                                    for i in range(1,32):
                                        t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^A[5]]<<i
                                    g0123456 = g012345 | (t<<32)
                                    if g0123456 not in fC_part[5]: continue
                            
                                    for A[6] in A_candidate[6]:
                                        t = table_g[A[6]]
                                        for i in range(1,64):
                                            t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^(((i>>5)&1)*A[5])^A[6]]<<i
                                        g01234567 = g0123456 | (t<<64)
                                        if g01234567 not in fC_part[6]: continue    
                                
                                        for A[7] in A_candidate[7]:
                                            t = table_g[A[7]]
                                            for i in range(1,128):
                                                t |= table_g[((i&1)*A[0])^(((i>>1)&1)*A[1])^(((i>>2)&1)*A[2])^(((i>>3)&1)*A[3])^(((i>>4)&1)*A[4])^(((i>>5)&1)*A[5])^(((i>>6)&1)*A[6])^A[7]]<<i
                                            g012345678 = g01234567 | (t<<128)
                                            if g012345678 not in fC_part[7]: continue
                                            if min_MC > MC:
                                                min_MC = MC
                                            Abox_set_g = affine_func[g]
                                            for Abox_g in Abox_set_g:
                                                Abox_f = []
                                                for i in range(len(Abox_g)-1):
                                                    h = Abox_g[i]
                                                    table_h = [(h>>i)&1 for i in range(1<<n)]
                                                    hA = 0
                                                    for i in range(1<<n):
                                                        xA = 0
                                                        for j in range(n):
                                                            xA ^= ((i>>j)&1)*A[j]
                                                        hA |= table_h[xA]<<i
                                                    Abox_f.append(hA)
                                                Abox_f.append(f)
                                                if f not in targets_Aboxes:
                                                    targets_Aboxes[f] = [Abox_f]
                                                else:
                                                    targets_Aboxes[f].append(Abox_f)
                                            f_flag = True
                                            break
                                        if f_flag: break
                                    if f_flag: break
                                if f_flag: break
                            if f_flag: break
                        if f_flag: break
                    if f_flag: break
                if f_flag: break
            if f_flag: break
    return min_MC,targets_Aboxes

def Bool_to_int(n,line_f):
    La_x = [0 for La in range(1<<n)]
    for La in range(1<<n):
        for x in range(1<<n):
            La_x[La] |= (bin(La&x).count('1')%2)<<x
    X = [La_x[1<<i] for i in range(n)]

    f = 0
    monos = line_f.split('+')
    for mono in monos:
        vars = list(map(int,re.findall(r'\d+',mono)))
        fmono = (1<<(1<<n))-1
        for var in vars:
            fmono &= X[var-1]
        f ^= fmono
    return f

def Read_Data_Abox_set(n,f):
    files = os.listdir('./Bool__')
    lines = []
    for filename in files:
        nums = list(map(lambda x:int(x),re.findall(r'\d+',filename)))
        if len(nums) < 4: continue
        # print(filename)
        if n == nums[0]:
            if f in nums:
                with open('./Bool__/'+filename,'r') as file:
                    lines = file.readlines()
                break
    # print(f,hex(f),filename)
    Abox_set = []
    for line in lines:
        Abox = list(map(lambda x:int(x),re.findall(r'\d+',line)))
        Abox_set.append(Abox)
    return Abox_set


def Read_Data_set_1(n,mc_dim):
    funcs_no_MC_dim = dict()
    for mc,dim in mc_dim:
        if dim > n: continue
        with open(f'./mc_dim/mc{mc}_dim{dim}.txt','r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            f = Bool_to_int(n,line)
            funcs_no_MC_dim[f] = [i,mc,dim]
    return funcs_no_MC_dim
