#P346 Computational Physics Lab
#My library 
#Vedansh Desai #2311198
import numpy as np
class MyComplex ( ) :
def init (self,real , imag = 0.0 ) :
self.r=real
self.i=imag
def displaycmplx ( self ) :
print( self.r ,”, ” , self . i , ” j ” , s e p=” ” )
def add cmplx ( s e l f , c1 , c 2 ) :
s e l f . r=c 1 . r+c 2 . r
s e l f . i=c 1 . i+c 2 . i
return MyComplex ( s e l f )
def s u b c m p l x ( s e l f , c1 , c 2 ) :
s e l f . r=c 1 . r−c 2 . r
s e l f . i=c 1 . i −c 2 . i
return MyComplex ( s e l f )
def mul cmplx ( s e l f , c1 , c 2 ) :
s e l f . r=c 1 . r ∗c 2 . r − c 1 . r ∗c 2 . i
s e l f . i=c 1 . i ∗c 2 . r + c 1 . r ∗c 2 . i
return MyComplex ( s e l f )
def mod cmplx ( s e l f ) :
return np . s q r t ( s e l f . r ∗∗2+ s e l f . i ∗∗2 )
2
