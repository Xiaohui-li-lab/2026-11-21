1 # !/ usr / bin / env python3
2 " " "
3 ligo_echo_analysis . py
4 Matched - filter upper - limit estimate for post - merger gravitational - wave
echoes .
5 " " "
6 import numpy as np
7 import matplotlib . pyplot as plt
8
9 def echo_template (t , A_echo , dt_echo , tau , f0 , n_echoes =5) :
10 h = np . zeros_like ( t )
11 for n in range ( n_echoes ) :
12 tn = ( n +1) * dt_echo
13 if tn > t [ -1]:
14 break
15 env = np . exp ( -(t - tn ) / tau ) * ( t >= tn )
16 h += A_echo * ( -1) ** n * env * np . sin (2* np . pi * f0 *( t - tn ) )
17 return h
18
19 def matched_filter_snr (h , template , psd , dt ) :
20 hf = np . fft . rfft ( h )
21 tf = np . fft . rfft ( template )
22 freq = np . fft . rfftfreq ( len ( h ) , dt )
23 psd_interp = np . interp ( freq , psd [: ,0] , psd [: ,1] , left = psd [0 ,1] , right =
psd [ -1 ,1])
24 df = freq [1] - freq [0]
25 snr_sq = 4* df * np . sum (( hf * tf . conj () ) . real / psd_interp )
26 norm = 4* df * np . sum ( np . abs ( tf ) **2 / psd_interp )
27 return snr_sq / np . sqrt ( norm )
28
9
29 def simulate_ligo_noise (N , dt , psd ) :
30 freq = np . fft . rfftfreq (N , dt )
31 psd_interp = np . interp ( freq , psd [: ,0] , psd [: ,1] , left = psd [0 ,1] , right =
psd [ -1 ,1])
32 amp = np . sqrt (0.5* psd_interp / dt )
33 phase = np . random . rand ( len ( freq ) ) *2* np . pi
34 return np . fft . irfft ( amp * np . exp (1 j * phase ) , n = N )
35
36 def run_echo_search ( psd_file = " aligo_psd . txt " , dt =1/4096 , T =8.0) :
37 psd = np . loadtxt ( psd_file )
38 N = int ( T / dt )
39 t = np . arange ( N ) * dt
40 noise = simulate_ligo_noise (N , dt , psd )
41 f0 = 250.0
42 tau = 0.1
43 dt_echo = 0.014
44 amplitudes = np . logspace ( -4 , -1 , 30)
45 snrs = []
46 for A in amplitudes :
47 tmpl = echo_template (t , A , dt_echo , tau , f0 )
48 snr = matched_filter_snr ( noise , tmpl , psd , dt )
49 snrs . append ( snr )
50 snrs = np . array ( snrs )
51 threshold = 5.0
52 A_95 = amplitudes [ np . argmax ( np . abs ( snrs ) > threshold ) ] if np . any ( np . abs (
snrs ) > threshold ) else amplitudes [ -1]
53 plt . figure ( figsize =(6 ,4) )
54 plt . semilogy ( amplitudes , np . abs ( snrs ) , label = ’| SNR | ’)
55 plt . axhline ( threshold , color = ’r ’ , ls = ’ -- ’ , label = ’ SNR =5 threshold ’)
56 plt . axvline ( A_95 , color = ’g ’ , ls = ’ -- ’ , label = f ’ 95% UL A ={ A_95 :.2 e } ’)
57 plt . xlabel ( ’ Echo amplitude $A_ {\\ rm echo } $ ’)
58 plt . ylabel ( ’ Matched - filter | SNR | ’)
59 plt . title ( ’ LIGO echo search upper limit ( toy simulation ) ’)
60 plt . legend ()
61 plt . tight_layout ()
62 plt . savefig ( ’ ligo_echo_constraint . png ’ , dpi =300)
63 print ( f " LIGO echo 95% upper limit : A_echo < { A_95 :.2 e } " )
64
65 if __name__ == " __main__ " :
66 run_echo_search ()
