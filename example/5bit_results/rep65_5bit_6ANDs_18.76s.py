# [1441853610, 11184640, 13421568, 1716013056, 15790080, 3223471296]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*6; y = [0]*5

	g[0] = (x[0]^x[2]^x[3]) & (x[0]^x[2]^x[4]) ^ x[2]
	g[1] = (x[0]) & (x[3]^x[4]) 
	g[2] = (x[1]) & (x[3]^x[4]) 
	g[3] = (x[0]^g[0]) & (x[0]^x[1]^g[1]) 
	g[4] = (x[2]) & (x[3]^x[4]) 
	g[5] = (x[0]^x[2]^x[3]^g[0]) & (x[1]^g[1]) 

	y[0] = g[0] 
	y[1] = x[0] ^ x[1] ^ x[3] ^ g[0] ^ g[1] ^ g[2] ^ g[3] ^ g[4] ^ g[5] 
	y[2] = x[2] ^ x[3] ^ x[4] ^ g[1] ^ g[3] ^ g[4] ^ g[5] 
	y[3] = x[0] ^ x[4] ^ g[0] ^ g[2] ^ g[3] ^ g[4] ^ g[5] 
	y[4] = x[0] ^ g[0] ^ g[1] ^ g[2] ^ g[4] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([0, 1, 2, 3, 4, 5, 8, 9, 6, 10, 16, 28, 7, 11, 31, 19, 12, 14, 20, 22, 13, 15, 27, 25, 17, 30, 29, 18, 21, 26, 23, 24])
if tuple(S)==tuple([0, 1, 2, 3, 4, 5, 8, 9, 6, 10, 16, 28, 7, 11, 31, 19, 12, 14, 20, 22, 13, 15, 27, 25, 17, 30, 29, 18, 21, 26, 23, 24]):
	print(True)
else:
	print(False)
