#################Point Class#################
class Point:
    def __init__(self, X, Y, F):
        ##X and Y are Finite Field Elements
        ##element is a binary list with least significant coef first
        ##ie x^5 + x^3 + x + 1 = [1,1,0,1,0,1]
        self.X = X.copy()
        self.Y = Y.copy()
        ##F is a polynomial that is irreducible over the Finite Field
        self.F = F.copy()
        ##Field is GF(2^n)
        self.n = len(self.F) - 1

    def getX(self):
        return self.X.copy()

    def getY(self):
        return self.Y.copy()

    def getF(self):
        ##get irreducible polynomial
        return self.F.copy()

    def getN(self):
        ##get n where GF(2^n)
        return self.n

    def padElements(self):
        ##Pad with 0s so elements are of same length
        l = self.n
        self.X += [0]*(l-len(self.X))
        self.Y += [0]*(l-len(self.Y))

    def isEqual(self, P):
        X1 = self.X.copy()
        Y1 = self.Y.copy()
        X2 = P.X.copy()
        Y2 = P.Y.copy()
        ##pad with 0s
        l = max(len(X1),len(X2),len(Y1),len(Y2))
        X1 += [0]*(l-len(X1))
        Y1 += [0]*(l-len(Y1))
        X2 += [0]*(l-len(X2))
        Y2 += [0]*(l-len(Y2))
        ##check if equal
        if (X1 == X2 and Y1 == Y2):
            return True
        else:
            return False

    def onCurve(self):
        return onCurve(self.getX(),self.getY(),self.getF())

    def add(self, P):
        return addPoints(self, P)

    def mult(self, k):
        return scalarMultPoint(self, k)

    def out(self):
        ##Returns point as binary string for printing
        x = ""
        y = ""
        for i in reversed(self.X.copy()):
            x += str(i)
        for i in reversed(self.Y.copy()):
            y += str(i)
        return("("+x+", "+y+")")
        
    def decOut(self):
        ##Returns point as decimal string for printing
        x = 0
        y = 0
        for i in range(len(self.X)):
            x += self.X[i] * 2**i
        for i in range(len(self.Y)):
            y += self.Y[i] * 2**i
        return("("+str(x)+", "+str(y)+")")

    def printPoly(self):
        poly = "("
        poly += printPoly(self.getX())
        poly += ", "
        poly += printPoly(self.getY())
        poly += ")"
        return poly

    def copy(self):
        return Point(self.X.copy(),self.Y.copy(),self.F.copy())

#################Finite Field Element Functions#################
def multFFE(A,B,F):
    ##Multiply two Finite Field Elements
    ##Output = A(x)B(x) mod F(x) mod 2
    ##Least significant coef first in array
    A = A.copy()
    B = B.copy()
    F = F.copy()
    if (F[-1] != 1):
        print("Check to ensure the degree of F(x) is correct")
    if (len(A)>len(F) or len(B)>len(F)):
        print("Check to ensure A(x) and B(x) have already been reduced")
    ##pad with 0s
    l = len(F)
    A += [0]*(l-len(A))
    B += [0]*(l-len(B))
    C = [0]*(l*2)
    ##Multiply (shift and add)
    for i in range(l):
        C = addFFE(C,B[i]*A)
        A = [0] + A
    ##mod F
    ##easier to do if most significant coef first
    C.reverse()
    F.reverse()
    ##check if C is 0
    if C.count(1) == 0:
        return [0]*(l-1)
    ##ensure they start with 1
    C = C[C.index(1):]
    F = F[C.index(1):]
    ##start dividing
    while(len(C) >= len(F)):
        C = addFFE(C,F)
        C = C[C.index(1):]
    C.reverse()
    return C

def inverseFFE(A,F):
    ##Inverse of a Finite Field Element
    ##(2^n)-2 = 2^(n-1) + 2^(n-2) + ... + 2
    m = len(F) - 1
    temp = multFFE(A,A,F)
    output = temp.copy()
    while (m > 2):
        temp = multFFE(temp,temp,F)
        output = multFFE(output,temp,F)
        m -= 1
    return output

def addFFE(A,B):
    ##Add two finite field elements (All this does is XOR them)
    ##Output = A(x) + B(x) mod 2
    ##Least significant coef first in array
    A = A.copy()
    B = B.copy()
    l = max(len(A),len(B))
    ##pad with 0s
    A += [0]*(l-len(A))
    B += [0]*(l-len(B))
    C = []
    for i in range(l):
        ## ^ is the XOR operator
        C.append(A[i]^B[i])
    return C

#################Point Functions#################
def addPoints(P1,P2):
    ##add points over elliptic curve y^2 + xy = x^3 + x^2 + 1
    ##where curve is over GF(2^n)
    F = P1.getF()
    if (P1.isEqual(P2)):
        ##calculate lambda
        t1 = multFFE(P1.getY(),inverseFFE(P1.getX(),F),F)
        t2 = P1.getX()
        s = addFFE(t1, t2)
        ##calculate X3
        t1 = multFFE(s,s,F)
        t2 = s
        t3 = [1]
        X3 = addFFE(t1,t2)
        X3 = addFFE(X3,t3)
    else:
        ##calculate lambda
        num = addFFE(P1.getY(),P2.getY())
        den = addFFE(P1.getX(),P2.getX())
        s = multFFE(num,inverseFFE(den,F),F)
        ##calculate X3
        t1 = addFFE(multFFE(s,s,F),s)
        t2 = addFFE(P1.getX(),P2.getX())
        t3 = [1]
        X3 = addFFE(t1,t2)
        X3 = addFFE(X3,t3)
    ##calculate Y3
    t1 = multFFE(addFFE(P1.getX(),X3),s,F)
    t2 = X3
    t3 = P1.getY()
    Y3 = addFFE(t1,t2)
    Y3 = addFFE(Y3,t3)
    ##create Point
    P3 = Point(X3,Y3,F)
    return P3

def scalarMultPoint(P,k):
    ##double and add method
    temp = P.copy()
    output = 0
    K = int2bin(k)
    for i in K:
        if (i == 1):
            if (output == 0):
                output = temp.copy()
            else:
                output = addPoints(output,temp)
        temp = addPoints(temp,temp)
    return output

def onCurve(X,Y,F):
    ##elliptic curve is of form y^2 + xy = x^3 + x^2 + 1 over GF(2^n)
    l = len(F) - 1
    left = addFFE(multFFE(Y,Y,F),multFFE(X,Y,F))
    left += [0]*(l-len(left))
    t2 = multFFE(X,X,F)
    t1 = multFFE(X,t2,F)
    t3 = [1]
    right = addFFE(t1,t2)
    right = addFFE(right,t3)
    right += [0]*(l-len(right))
    if (left == right):
        return True
    else:
        return False

#################Elliptic Curve Functions#################
def orderPoint(A,P,N):
    ##Find the order of point A on EC
    ##N is the number of points on EC
    ##P is a primitive point on the EC
    ##elliptic curve is of form y^2 + xy = x^3 + x^2 + 1 over GF(2^n)
    temp = P.copy()
    ##Find k where kP = a
    for k in range(1,N):
        if temp.isEqual(A):
            break
        temp = addPoints(temp,P)
    ##Find order
    for i in range(1,N+1):
        if (k*i%N==0):
            return i
    return 0

def orderSimple(i,E):
    for x in range(1,E+1):
        if(i*x%E==0):
            return x
    return 0

def numberPoints(F):
    ##get the number of points on the EC
    ##F is the irreducible polynomial
    ##elliptic curve is of form y^2 + xy = x^3 + x^2 + 1 over GF(2^n)
    pointList = []
    l = len(F) - 1
    i = 0
    for xi in range(2**l):
        X = int2bin(xi)
        for yi in range(2**l):
            Y = int2bin(yi)
            ##Test all points to check if on curve
            if (onCurve(X,Y,F)):
                i += 1
                pointList.append(Point(X,Y,F))
        if (xi%30 == 0):
            ##Print Progress
            print("Found",xi,"points so far...")
    ##Don't forget Point O
    i += 1
    print("Number of Points =",str(i))
    return i

def maximalPoints(P,N):
    ##Primitive point P
    ##prime p
    ##value a
    ##N = number of points on curve E
    P.padElements()
    testPoint = P.copy()
    i = 1
    while(True):
        order = orderSimple(i,N)
        ##order = orderPoint(testPoint,P,N)
        if (testPoint.onCurve() == False):
            print("Error: Point",testPoint.out(),"not on curve!")
            break
        if (order == N):
            ##Print Maximal Points
            print (i,"P =",testPoint.out(),"ord =",str(order))
        if(testPoint.getX() == P.getX() and i > 1 and i < N-1):
            ##The program should end with the correct number of points
            ##where the last point is the negative of the first point
            ##since adding the two results in Point 0
            print("Error: did not start with primitive element")
            break
        if(i==N-1):
            break
        i += 1
        testPoint = addPoints(testPoint,P)
        testPoint.padElements()
    return 0

#################List Functions#################
##These Functions are mostly for making the inputs/outputs more readable
##helps prevents careless mistakes when checking my work
def printPoly(P):
    ##print the polynomial in a readable fashion
    ##Least significant coef first in array
    output = ""
    if (P[0] == 1):
        output = "1"
    for i in range (1,len(P)):
        if (P[i] == 1):
            output = "X^" + str(i) + " + " + output
    if (output[-2] == "+"):
        output = output[0:-3]
    return output

def exp2bin(exp):
    ##takes list of exponents and converts to binary list
    ##ie x^5 + x^3 + x + 1 == [5,3,1,0] => [1,1,0,1,0,1]
    output = [0] * (max(exp) + 1)
    for i in exp:
        output[i] = 1
    return output

def str2bin(s):
    ##Input: "10111" Output: [1,1,1,0,1]
    s = list(s)
    s.reverse()
    for i in range(len(s)):
        s[i] = int(s[i])
    return s

def hex2bin(h):
    ##takes hex string and converts to binary list in reverse order
    ##ie "50" => 0101 0000 => [0,0,0,0,1,0,1,0]
    output = []
    hexTable = ['0000', '0001', '0010', '0011',
                '0100', '0101', '0110', '0111',
                '1000', '1001', '1010', '1011',
                '1100', '1101', '1110', '1111']
    for c in h:
        output += list(hexTable[int(c,16)])
    output.reverse()
    for i in range(len(output)):
        output[i] = int(output[i])
    return output

def int2bin(i):
    ##takes a decimal integer and converts to a binary list in reverse order
    if (i == 0):
        return [0]
    ##convert to hex string first
    h = hex(i)[2:]
    output = hex2bin(h)
    ##remove trailing 0s
    output.reverse()
    output = output[output.index(1):]
    output.reverse()
    return output

#################ECDHKE#################
##irreducible polynomial
F = exp2bin([9,8,0])
##Primitive Point
Xp = exp2bin([1,0])
Yp = exp2bin([5,4,3])
P = Point(Xp,Yp,F)
P.padElements()
##number of points on EC
##Function checks every point and takes about a minute to run
##Uncomment to run it
##N = numberPoints(F)
N = 518
##Alice and Bob's private key
a = 13
b = 31
print("Primitive point P =",P.out())
##Step 1
print("Step 1")
A = scalarMultPoint(P,a)
B = scalarMultPoint(P,b)
print("Alice computes 13P = A =",A.out())
print("And sends point A to Bob")
print("Bob computes 31P = B =",B.out())
print("And sends point B to Alice")
##Step 2
print("Step 2")
A2 = scalarMultPoint(B,a)
B2 = scalarMultPoint(A,b)
print("Alice receives B =",B.out(),"and computes 13B =",A2.out())
print("Bob receives A =",A.out(),"and computes 31A =",B2.out())
##Confirm
C = scalarMultPoint(P,a*b)
print("Both parties have the shared key",C.out())
##get Maximal Points
print("Maximal Points")
m = maximalPoints(P,N)

#################Unused Functions#################
##I wrote these functions using alternate algorithms
##They were created to help verify the correctness of my outputs
##They are not used in the program
def multFFE2(A,B,F):
    ##algorithm from slides
    ##Output = A(x)B(x) mod F(x) mod 2
    ##Least significant coef first in array
    A = A.copy()
    B = B.copy()
    F = F.copy()
    if (F[-1] != 1):
        print("Check to ensure the degree of F(x) is correct")
    if (len(A)>len(F) or len(B)>len(F)):
        print("Check to ensure A(x) and B(x) have already been reduced")
    ##pad with 0s
    l = len(F)
    A += [0]*(l-len(A))
    B += [0]*(l-len(B))
    ##irreducible polynomial replacement
    f = F[0:-1]
    ##begin multiplication
    C = [[]]*l
    ##First Step
    C[0] = [0]+B[-1]*A
    C[0] += [0]*(l-len(C[0]))
    if (C[0][l-1] == 1):
        C[0] = addFFE(C[0],f)
    C[0] = C[0][0:l-1]
    ##iterate
    for i in range(1,l-1):
        C[i] = [0]+addFFE(C[i-1],B[-(i+1)]*A)
        C[i] += [0]*(l-len(C[i]))
        if (C[i][l-1] == 1):
            C[i] = addFFE(C[i],f)
        C[i] = C[i][0:l-1]
    ##final step
    C[l-1] = addFFE(C[l-2],B[0]*A)
    C[l-1] += [0]*(l-len(C[l-1]))
    if (C[l-1][l-1] == 1):
        C[l-1] = addFFE(C[l-1],f)
    C[l-1] = C[l-1][0:l-1]
    return C[l-1]

def inverseFFEBrute(A,F,m):
    ##Output is inverse of A(x) mod F(x) in GF(2^m)
    A = A.copy()
    F = F.copy()
    e = 2**m-2
    temp = A.copy()
    for i in range(e-1):
        temp = multFFE(temp,A,F)
    if  (multFFE(A,temp,F) != [1]):
        print ("Inverse Error")
    return temp

def inverseFFERecursive(A,F,m):
    ##not needed, just wrote this method for fun
    ##(2^n)-2 = 2^(n-1) + 2^(n-2) + ... + 2
    A = A.copy()
    m -= 1
    if (m > 1):
        temp = A.copy()
        for i in range(m):
            temp = multFFE(temp,temp,F)
        return multFFE(temp,inverseRecursive(A,F,m),F)
    else:
        return multFFE(A,A,F)
        
def scalarMultPointSimple(P,k):
    output = P.copy()
    while(k>=2):
        output = addPoints(output,P)
        k -= 1
    return output
    
def powFFE(A,F,n):
    ##compute A^n
    output = A.copy()
    while (n >= 2):
        output = multFFE(output,A,F)
        n -= 1
    return output

def gcd (a,b):
    if a == 0 or b == 0:
        return 0
    while a != 0 and b != 0:
        if a > b:
            a %= b
        else:
            b %= a
    if a == 0:
        return b
    else:
        return a
