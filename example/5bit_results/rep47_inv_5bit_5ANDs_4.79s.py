# [2526412800, 4032482544, 15790080, 2768256000, 178261184]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*5; y = [0]*5

	g[0] = (x[0]^x[1]^x[2]) & (x[4]) 
	g[1] = (x[3]^x[4]) & (x[1]^x[2]^g[0]) ^ x[2]
	g[2] = (x[2]) & (x[3]^x[4]) 
	g[3] = (x[3]) & (x[1]^x[4]^g[0]^g[2]) 
	g[4] = (x[2]^x[3]) & (x[1]^g[0]^g[2]) 

	y[0] = x[0] ^ g[0] ^ g[4] 
	y[1] = x[1] ^ x[2] ^ g[0] ^ g[1] ^ g[2] ^ g[3] 
	y[2] = g[1] 
	y[3] = x[2] ^ x[3] ^ g[1] ^ g[3] 
	y[4] = x[4] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 15, 14, 10, 11, 12, 13, 16, 28, 17, 29, 21, 24, 20, 25, 18, 27, 19, 26, 31, 23, 30, 22])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 15, 14, 10, 11, 12, 13, 16, 28, 17, 29, 21, 24, 20, 25, 18, 27, 19, 26, 31, 23, 30, 22]):
	print(True)
else:
	print(False)
