# [2694881440, 969340170, 574899234, 806105100, 94394880, 89169930, 1145317416]
S = []
for X in range(32):
	x=[(X>>i)&1 for i in range(5)]
	g = [0]*7; y = [0]*5

	g[0] = (x[0]) & (x[2]) 
	g[1] = (x[1]) & (x[4]) ^ x[0]^ x[3]^ g[0]
	g[2] = (x[0]^x[1]) & (x[0]^x[3]^x[4]) 
	g[3] = (x[1]^x[2]) & (x[1]^x[3]) 
	g[4] = (x[0]^x[2]^x[4]) & (x[3]^g[0]) 
	g[5] = (x[0]^x[2]^x[3]) & (x[0]^x[4]) 
	g[6] = (x[0]^x[4]) & (x[1]^g[0]) 

	y[0] = g[1] ^ 1 
	y[1] = x[3] ^ x[4] ^ g[0] ^ g[1] ^ g[4] ^ g[5] ^ g[6] 
	y[2] = x[0] ^ x[1] ^ x[3] ^ x[4] ^ g[2] ^ g[3] ^ g[4] ^ g[6] 
	y[3] = x[4] ^ g[0] ^ g[3] 
	y[4] = x[0] ^ x[2] ^ g[0] ^ g[1] ^ g[2] ^ g[3] ^ g[4] ^ g[5] ^ g[6] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4))

print(S)
print([1, 0, 25, 26, 17, 29, 21, 27, 20, 5, 4, 23, 14, 18, 2, 28, 15, 8, 6, 3, 13, 7, 24, 16, 30, 9, 31, 10, 22, 12, 11, 19])
if tuple(S)==tuple([1, 0, 25, 26, 17, 29, 21, 27, 20, 5, 4, 23, 14, 18, 2, 28, 15, 8, 6, 3, 13, 7, 24, 16, 30, 9, 31, 10, 22, 12, 11, 19]):
	print(True)
else:
	print(False)
