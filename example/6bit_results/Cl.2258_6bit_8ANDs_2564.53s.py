# [57645195634461696, 14236000717916561610, 10377992699400720390, 11914835766936207360, 14051231095106780160, 14051231095097967360, 11551909311015833760, 11024980114793326080]
S = []
for X in range(64):
	x=[(X>>i)&1 for i in range(6)]
	g = [0]*8; y = [0]*6

	g[0] = (x[1]) & (x[3]^x[4]) 
	g[1] = (x[2]^x[3]^x[4]) & (x[0]^x[1]^x[5]) ^ x[0]^ x[4]^ g[0]
	g[2] = (x[0]^x[1]^x[3]) & (x[0]^x[1]^x[2]^x[4]^x[5]) 
	g[3] = (x[0]^x[2]^x[3]) & (x[4]) 
	g[4] = (x[1]^x[2]^x[4]) & (x[1]^x[2]^x[3]^x[5]) ^ x[1]^ x[2]^ x[4]
	g[5] = (x[1]^x[2]^x[3]) & (x[1]^x[2]^x[4]^x[5]) ^ x[1]^ x[2]^ x[3]
	g[6] = (x[2]) & (x[0]^x[3]^x[5]) 
	g[7] = (x[3]) & (x[0]^x[1]^x[5]) 

	y[0] = g[1] 
	y[1] = x[3] ^ x[5] ^ g[1] ^ g[2] ^ g[4] ^ g[5] ^ g[6] 
	y[2] = x[2] ^ x[4] ^ g[4] ^ g[7] 
	y[3] = x[4] ^ x[5] ^ g[3] ^ g[5] 
	y[4] = x[1] ^ x[3] ^ x[4] ^ g[0] ^ g[1] ^ g[2] ^ g[3] ^ g[5] 
	y[5] = x[3] ^ g[0] ^ g[5] 
	S.append((y[0]<<0)|(y[1]<<1)|(y[2]<<2)|(y[3]<<3)|(y[4]<<4)|(y[5]<<5))

print(S)
print([0, 1, 2, 3, 4, 6, 7, 5, 8, 12, 16, 20, 32, 39, 57, 62, 9, 17, 21, 13, 41, 50, 52, 47, 55, 42, 49, 44, 59, 37, 60, 34, 10, 25, 38, 53, 35, 51, 14, 30, 61, 43, 11, 29, 56, 45, 15, 26, 22, 28, 36, 46, 27, 18, 40, 33, 23, 24, 63, 48, 54, 58, 31, 19])
if tuple(S)==tuple([0, 1, 2, 3, 4, 6, 7, 5, 8, 12, 16, 20, 32, 39, 57, 62, 9, 17, 21, 13, 41, 50, 52, 47, 55, 42, 49, 44, 59, 37, 60, 34, 10, 25, 38, 53, 35, 51, 14, 30, 61, 43, 11, 29, 56, 45, 15, 26, 22, 28, 36, 46, 27, 18, 40, 33, 23, 24, 63, 48, 54, 58, 31, 19]):
	print(True)
else:
	print(False)
