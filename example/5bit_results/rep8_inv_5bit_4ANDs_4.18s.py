# [1718004394, 2863267840, 4042260480, 4278190080]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*4; y = [0]*5

	g[0] = (x[1]) & (x[4]) ^ x[0]
	g[1] = (x[0]) & (x[4]) 
	g[2] = (x[2]) & (x[4]) 
	g[3] = (x[3]) & (x[4]) 

	y[0] = g[0] 
	y[1] = x[0] ^ x[1] ^ g[0] ^ g[1] 
	y[2] = x[2] ^ g[3] 
	y[3] = x[3] ^ g[2] ^ g[3] 
	y[4] = x[4] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 17, 18, 28, 31, 29, 30, 20, 23, 21, 22, 24, 27, 25, 26])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 17, 18, 28, 31, 29, 30, 20, 23, 21, 22, 24, 27, 25, 26]):
	print(True)
else:
	print(False)
