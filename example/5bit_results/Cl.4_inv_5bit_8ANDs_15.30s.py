# [3221425200, 1352683680, 3420464, 898825552, 21266710, 2190284240, 2351042944, 42341958]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*8; y = [0]*5

	g[0] = (x[1]^x[3]) & (x[1]^x[2]^x[4]) ^ x[1]^ x[2]^ x[4]
	g[1] = (x[2]) & (x[0]^x[3]) 
	g[2] = (x[3]^x[4]^g[0]) & (x[0]^x[1]^x[2]^x[3]^g[1]) 
	g[3] = (x[0]^g[1]) & (x[0]^x[4]^g[2]) ^ x[0]^ x[2]^ x[4]^ g[0]^ g[2]
	g[4] = (x[0]^x[1]^x[3]^g[0]) & (x[0]^x[1]^g[1]^g[3]) 
	g[5] = (x[0]^x[2]^g[0]^g[1]^g[2]) & (x[0]^x[1]^x[2]^x[3]^x[4]^g[3]^g[4]) 
	g[6] = (x[1]^x[2]^g[0]^g[1]^g[2]) & (x[2]^x[4]^g[1]^g[5]) 
	g[7] = (x[0]^x[1]^g[1]) & (x[2]^x[3]^g[4]^g[6]) 

	y[0] = x[1] ^ x[3] ^ x[4] ^ g[1] ^ g[2] ^ g[3] ^ g[4] ^ g[5] ^ g[6] 
	y[1] = x[1] ^ x[2] ^ x[4] ^ g[0] ^ g[3] ^ g[5] ^ g[6] 
	y[2] = g[3] 
	y[3] = x[0] ^ x[1] ^ x[3] ^ g[0] ^ g[2] ^ g[7] 
	y[4] = x[2] ^ x[3] ^ x[4] ^ g[0] ^ g[1] ^ g[2] ^ g[4] ^ g[6] ^ g[7] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 8, 5, 16, 6, 9, 20, 7, 12, 29, 22, 30, 10, 13, 11, 19, 31, 18, 27, 15, 28, 25, 21, 26, 14, 23, 24, 17])
if tuple(S)==tuple([0, 1, 2, 3, 4, 8, 5, 16, 6, 9, 20, 7, 12, 29, 22, 30, 10, 13, 11, 19, 31, 18, 27, 15, 28, 25, 21, 26, 14, 23, 24, 17]):
	print(True)
else:
	print(False)
