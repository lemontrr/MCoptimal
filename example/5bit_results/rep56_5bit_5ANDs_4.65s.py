# [4042312908, 2852170240, 3233857728, 3435921408, 1437204480]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*5; y = [0]*5

	g[0] = (x[1]^x[2]) & (x[4]) ^ x[1]
	g[1] = (x[0]) & (x[3]) 
	g[2] = (x[1]) & (x[2]) 
	g[3] = (x[1]) & (x[4]) 
	g[4] = (x[0]^x[3]) & (x[4]) 

	y[0] = x[0] ^ x[1] ^ x[4] ^ g[0] ^ g[1] ^ g[2] ^ g[3] 
	y[1] = g[0] 
	y[2] = x[2] ^ g[3] 
	y[3] = x[3] ^ x[4] ^ g[1] 
	y[4] = g[1] ^ g[4] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 7, 6, 8, 16, 10, 18, 12, 20, 15, 23, 9, 24, 13, 28, 14, 31, 11, 26, 17, 25, 21, 29, 22, 30, 19, 27])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 7, 6, 8, 16, 10, 18, 12, 20, 15, 23, 9, 24, 13, 28, 14, 31, 11, 26, 17, 25, 21, 29, 22, 30, 19, 27]):
	print(True)
else:
	print(False)
