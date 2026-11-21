# !/ usr / bin / env python3
2 " " "
3 point_decomposition . py
4 Reproducible toy model of the consciousness / composite - point subdivision .
5 " " "
6 import numpy as np
7 import pandas as pd
8
9 def decompose_point ( epsilon0 =1.0 , levels =11 , alpha =0.5) :
10 rows = []
11 for n in range ( levels ) :
12 epsilon = epsilon0 * alpha ** n
13 nodes = 2** n
14 dim = 2** n
4
15 entropy = np . log ( dim ) if dim > 0 else 0.0
16 rows . append ({
17 ’ level ’: n ,
18 ’ scale ’: epsilon ,
19 ’ nodes ’: nodes ,
20 ’ hilbert_dim ’: dim ,
21 ’ entropy ’: entropy
22 })
23 return pd . DataFrame ( rows )
24
25 if __name__ == ’ __main__ ’:
26 df = decompose_point ( epsilon0 =1.0 , levels =11 , alpha =0.5)
27 df . to_csv ( ’ point_decomposition . csv ’ , index = False )
28 print ( df . to_string ( index = False ) )
29 slope , intercept = np . polyfit ( np . log (1/ df [ ’ scale ’ ]) , df [ ’ entropy ’] , 1)
30 expected = np . log (2) / np . log (1/0.5)
31 print ( f ’\ nFitted slope = { slope :.4 f } ( expected = { expected :.4 f }) ’)
