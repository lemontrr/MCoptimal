# [3275504700, 3224383536, 2583718400, 2590638080, 1315900]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*5; y = [0]*5

	g[0] = (x[3]) & (x[4]) ^ x[1]^ x[2]
	g[1] = (x[2]) & (g[0]) 
	g[2] = (x[3]) & (x[0]^x[2]^g[1]) 
	g[3] = (x[4]) & (x[0]^x[2]^g[1]) 
	g[4] = (x[1]^x[2]) & (g[0]^g[2]^g[3]) 

	y[0] = x[0] ^ x[2] ^ x[4] ^ g[1] ^ g[2] 
	y[1] = x[2] ^ g[4] 
	y[2] = x[2] ^ g[0] ^ g[4] 
	y[3] = x[3] ^ x[4] ^ g[2] 
	y[4] = x[1] ^ x[2] ^ g[0] ^ g[2] ^ g[3] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 7, 6, 8, 16, 10, 20, 12, 18, 22, 14, 9, 24, 11, 28, 13, 26, 30, 15, 21, 29, 17, 25, 31, 23, 19, 27])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 7, 6, 8, 16, 10, 20, 12, 18, 22, 14, 9, 24, 11, 28, 13, 26, 30, 15, 21, 29, 17, 25, 31, 23, 19, 27]):
	print(True)
else:
	print(False)
