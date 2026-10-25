# SIPI 學期總回顧

依照[原始重點回顧問題](instructions/重點回顧問題.txt)的順序，逐題回答 **89 題、13 個主題**。保留英文原題，解答使用繁體中文，涵蓋傅立葉轉換、數位訊號頻譜、眼圖、等化器、差動傳輸、串音、抖動、PAM4，以及展頻技術。

各題末的 `[Lxx, p.N]` 對應文末[講義對照](#lecture-sources)；頁碼採 PDF 從封面起算的頁次。公式中的假設與適用條件列在相關解答內。這份回顧可搭配[各次作業](../../README.md)的程式與模擬結果閱讀。

## 閱讀順序

| 題號 | 主題 | 對照講義 |
| --- | --- | --- |
| [Q01–Q07](#topic-01) | 傅立葉轉換 | L01 |
| [Q08–Q12](#topic-02) | 反傅立葉轉換 | L02 |
| [Q13–Q15](#topic-03) | 分段卷積 | L03 |
| [Q16–Q22](#topic-04) | 數位訊號的頻譜 | L04 |
| [Q23–Q26](#topic-05) | 眼圖基礎 | L05 |
| [Q27–Q31](#topic-06) | 疊加、ISI 與 PDA | L06 |
| [Q32–Q43](#topic-07) | CTLE、FFE 與 DFE | L07–L09 |
| [Q44–Q51](#topic-08) | 差動訊號一 | L10 |
| [Q52–Q59](#topic-09) | 串音 | L11 |
| [Q60–Q71](#topic-10) | 抖動與 BER | L12–L13 |
| [Q72–Q78](#topic-11) | 差動訊號二 | L14 |
| [Q79–Q82](#topic-12) | PAM4 | L15 |
| [Q83–Q89](#topic-13) | Scrambler 與 SSC | L16 |

## 符號約定

- 傅立葉轉換採 $X(f)=\int x(t)e^{-j2\pi ft}\,dt$；DFT 不含正規化因子，IDFT 含 $1/N$。$\operatorname{sinc}(u)=\sin(\pi u)/(\pi u)$。
- 除非另有說明，電壓 PSD 為雙邊頻譜，單位是 $\mathrm{V^2/Hz}$；跨電阻 $R$ 的功率 PSD 須再除以 $R$。單邊頻譜、峰值振幅、RMS 值不能直接互換。
- NRZ 的一個 UI 為位元時間 $T_b$；PAM4 的一個 UI 為符元時間 $T_s$，每個符元帶兩個位元。比較兩者時須先固定比較的是位元率、符元率或總擺幅。
- $S_{21}$ 是波的傳輸係數。在兩埠參考阻抗相同且源端、負載端均匹配時，若輸入定義為訊號源的開路電壓，輸出電壓轉移為 $S_{21}/2$；輸入已定義為入射電壓波時，不可再除以二。未匹配時還須計入反射。
- 差動通道的四埠例子採輸入對 $(1,3)$、輸出對 $(2,4)$；混模排序為 $[c_1,d_1,c_2,d_2]$。更換埠順序時，須同步改變轉換矩陣。

<a id="topic-01"></a>

## Fourier transform

<a id="q01"></a>

### Q01 What is the difference between DFT and FFT?

DFT（discrete Fourier transform）是**數學轉換的定義**。對 $N$ 個樣本 $x[n]$，本篇採用未正規化的正向轉換：$X[k]=\sum_{n=0}^{N-1}x[n]e^{-j2\pi kn/N}$。它把有限長度的時間樣本，轉成 $N$ 個等間距的頻率樣本；若把索引延伸到這個範圍以外，兩個序列都以 $N$ 為週期。[L01, p.2]

FFT（fast Fourier transform）是**計算同一個 DFT 的快速演算法**，並非另一種頻譜。直接依定義計算約需 $O(N^2)$ 次運算，常見 FFT 約需 $O(N\log N)$；在相同定義與正規化下，兩者結果只可能因浮點誤差略有差異。Radix-2 FFT 適合 $N=2^m$，但 FFT 並不限定只能處理二的冪次長度。[L01, p.2; L03, pp.2, 13]

<a id="q02"></a>

### Q02 What's the procedure to compute FT via DFT?

先取得等間距樣本 $x[n]=x(t_0+n\Delta t)$，共 $N_0$ 點，觀察長度定為 $T_0=N_0\Delta t$。取樣率 $f_s=1/\Delta t$ 必須足以涵蓋要分析的頻帶，且取樣前需避免頻帶外成分混疊。依動態範圍與頻率解析度需求選擇窗函數 $w[n]$，先對**原始的 $N_0$ 點**加窗，再視需要補零到 $M\geq N_0$ 點。對窗化、截斷訊號的連續 FT，可用下式近似：

$$
X_w(f_k)\approx \Delta t\,e^{-j2\pi f_k t_0}
\operatorname{FFT}_M\{w[n]x[n]\}[k],
\qquad \Delta f=\frac{1}{M\Delta t}.
$$

頻率軸必須依 FFT 索引標示；後半段索引代表負頻率，可用 `fftfreq`／`fftshift` 正確排列。若估計正弦波振幅，還要補償窗的 coherent gain，並區分單邊與雙邊表示；若估計隨機訊號 PSD，則使用 Q19 的功率正規化。加窗造成的頻譜平滑不能只靠除以一個增益完全還原，因此上式嚴格對應的是 $w(t)x(t)$ 的 FT。[L01, pp.8–10, 18–19]

<a id="q03"></a>

### Q03 What's the purpose of zero padding in DFT?

時域補零會讓 DFT 在**相同的 DTFT 曲線上取更多頻率樣本**。原本 $N_0$ 點的頻率格距為 $f_s/N_0$，補零到 $M$ 點後變為 $f_s/M$，因此曲線更平滑，較容易讀取峰值位置、零點與主瓣形狀，也方便不同資料使用共同的頻率格點。[L01, p.10]

補零沒有增加觀察到的資料，不能把原本無法分辨的兩個鄰近頻率真正分開。分辨能力仍由原始觀察時間 $T_0=N_0\Delta t$ 和窗的主瓣寬度決定，不能把補零後較小的 $\Delta f$ 當作較高的物理解析度。講義例子中，原始資料的解析度約為 $12.1\ \mathrm{MHz}$，即使補零把格距降到 $2\ \mathrm{MHz}$，解析度也不會變成 $2\ \mathrm{MHz}$。[L01, pp.10, 13]

<a id="q04"></a>

### Q04 What's the relation between dt and df in DFT?

若實際執行的是 $M$ 點 DFT，則

$$
\Delta t\,\Delta f=\frac{1}{M},\qquad
f_s=\frac{1}{\Delta t},\qquad
\Delta f=\frac{f_s}{M}.
$$

必須分清楚**原始資料長度 $N_0$**與**轉換長度 $M$**：未補零時 $M=N_0$，所以 $\Delta f=1/T_0$；補零後雖然 $\Delta f=1/(M\Delta t)$ 變小，真正觀察時間仍是 $T_0=N_0\Delta t$。$M\Delta t$ 是 DFT 所隱含的時間週期，包含人為加入的零值，不能當作新增的量測時間。[L01, pp.9–10]

<a id="q05"></a>

### Q05 What are the effects of time-domain truncation?

只保留長度 $T_0$ 的波形，等於把原訊號乘上一個矩形窗；在頻域則是原頻譜與矩形窗的 sinc 頻譜做卷積。例如 $0\leq t<T_0$ 的矩形窗具有 $W(f)=T_0\operatorname{sinc}(fT_0)e^{-j\pi fT_0}$，所以原本很窄的頻譜線會被展寬，鄰近成分被混合，稱為 **spectral blurring**；其旁瓣會把強訊號的能量散到其他頻率，稱為 **spectral leakage**，可能遮住很小的訊號。[L01, pp.14–15]

延長有效觀察時間能縮窄主瓣；若波形在截斷前已衰減到足夠小，截斷誤差也會較小。講義以矩形窗的 **6 dB 全主瓣寬度**定義解析度，得到約 $1.21/T_0$；這是特定窗與寬度準則，不能直接套用到所有窗。截斷洩漏與取樣不足造成的 aliasing 是不同問題，單純補零不能消除其中任何一種。[L01, pp.8, 15–17]

<a id="q06"></a>

### Q06 What's the purpose of using window function in DFT?

窗函數使有限記錄的兩端較平順地衰減，降低截斷邊界的不連續程度，因而壓低頻譜旁瓣，減少強成分洩漏到弱成分附近。頻域上，$\mathcal F\{w(t)x(t)\}=W(f)*X(f)$；選窗就是在選擇這個平滑核的主瓣與旁瓣特性。[L01, p.16]

代價通常是主瓣變寬，兩個鄰近頻率更難分開，且量到的振幅或雜訊功率必須做正確的增益補償。例如講義的 Hann 窗可明顯降低遠處旁瓣，但 6 dB 主瓣寬度約由矩形窗的 $1.21/T_0$ 增加到 $2/T_0$。因此選窗應由「要分開多近的頻率」與「要在多強的訊號旁看到多弱的訊號」共同決定。[L01, pp.16–18]

<a id="q07"></a>

### Q07 What is the gain of a window?

「窗的增益」必須說明是振幅還是功率。對長度 $N_0$ 的窗，**coherent gain** 為 $G_c=\frac{1}{N_0}\sum_n w[n]$，描述落在頻率格點上的正弦波，其 FFT 峰值振幅相對矩形窗縮小多少；估計正弦波振幅時以它補償。**Incoherent／power gain** 為 $G_p=\frac{1}{N_0}\sum_n w^2[n]$，描述白雜訊等隨機訊號的平均頻譜功率縮放；PSD 正規化應用 $G_p$，不能直接拿 $G_c^2$ 代替。若以 RMS 振幅增益表示，則為 $\sqrt{G_p}$；講義這裡的 incoherent gain 採功率定義。[L01, p.19; L04, pp.15–16]

兩者也決定窗的等效雜訊頻寬（ENBW）：

$$
B_{\mathrm{ENBW}}
=f_s\frac{\sum_n w^2[n]}{\left(\sum_n w[n]\right)^2}
=\frac{f_s}{N_0}\frac{G_p}{G_c^2}.
$$

例如 periodic Hann 窗有 $G_c=1/2$、$G_p=3/8$，因此 ENBW 為 $1.5f_s/N_0$；有限長度的 symmetric Hann 數值會略有差異，應直接由實際窗樣本計算。補零不改變這些加總，也不改變以 Hz 表示的 ENBW。ENBW、coherent gain 與 3 dB／6 dB 主瓣寬度分別描述不同特性。[L01, pp.17–19; L04, pp.15–16; W01]

<a id="topic-02"></a>

## Inverse Fourier transform

<a id="q08"></a>

### Q08 What's the spectrum symmetry property?

若時域訊號 $x(t)$ 為實數，其 FT 具有 Hermitian（共軛）對稱性：$X(-f)=X^*(f)$。因此實部是偶函數、虛部是奇函數、幅度為偶函數，而相位在模 $2\pi$ 意義下為奇對稱。注意這不是說 $X(-f)=X(f)$；若只鏡射幅度而不共軛相位，通常會重建出錯誤波形。[L02, p.3]

對實數序列的 $M$ 點 DFT，對應關係是 $X[M-k]=X^*[k]$，索引按 $M$ 取模。DC 點 $X[0]$ 必須為實數；若 $M$ 為偶數，Nyquist 點 $X[M/2]$ 也必須為實數，兩者都不能當作一對正、負頻率重複加入。一般複數時域訊號不必滿足這個對稱性。[L02, p.4]

<a id="q09"></a>

### Q09 What's the procedure of computing IFT via IDFT?

先確定複數頻譜的頻率間距 $\Delta f$ 均勻、包含 DC，且頻率上限 $f_{\max}$ 足以保留關心的波形細節。選定所需的 $\Delta t$，令 $M=1/(\Delta f\Delta t)$ 為可用的整數，並確保 $\Delta t\leq1/(2f_{\max})$。若資料缺少 DC，需依物理模型合理補值；直接填零可能改變最終電壓與低頻行為。重建波形以 $1/\Delta f$ 為週期，因此這個時間範圍要足以容納響應，避免週期回繞。[L02, pp.5–8]

一種做法是依 Q08 建立包含正、負頻率的完整 Hermitian 頻譜 $X_H[k]$，適當位置補零，再計算

$$
x[n]\approx M\Delta f\,\operatorname{IFFT}_M\{X_H[k]\},
\qquad t_n=n\Delta t.
$$

講義採用等價的「只存正頻率」方法：把 DC 除以 2，將資料補零到 $M$ 點後，取 $2M\Delta f\,\operatorname{Re}\{\operatorname{IFFT}(X_+)\}$；若包含恰好位於 Nyquist 的點，也要避免將它重複計數。這兩種方法的縮放不同，不能混用。需要降低頻寬截斷振鈴時，先對頻譜乘以合適的頻域窗；若使用完整 Hermitian 頻譜或 `irfft`，DC 本身不再先除以 2。[L02, pp.6–7, 12, 15]

<a id="q10"></a>

### Q10 What's the purpose of zero padding in IDFT?

固定 $\Delta f$，將轉換長度增加到 $M$，可把時間格距縮小為 $\Delta t=1/(M\Delta f)$，讓同一個帶限重建波形有更多時間樣本。這對讀取過零時間、上升時間、脈衝峰值或眼圖很有幫助，但時間週期 $M\Delta t=1/\Delta f$ 不變，原始頻率上限也沒有增加。[L02, p.7]

因此補零提供的是時間插值，無法還原沒有量到的高頻資訊，也不能提升由 $f_{\max}$ 與頻域窗決定的真實時間解析度。使用完整雙邊頻譜時，零必須放在正、負頻率之間的未量測高頻區，維持正確的頻率索引與共軛對稱；把已排好的雙邊頻譜直接在尾端補零，會改變原有頻率的物理位置。[L02, pp.7, 9, 12–13]

<a id="q11"></a>

### Q11 What are the effects of spectrum truncation?

把 $|f|>f_{\max}$ 的頻譜直接設為零，相當於乘上頻域矩形窗。若 $F=2f_{\max}$，其時域效應為

$$
x_F(t)=x(t)*\left[F\operatorname{sinc}(Ft)\right].
$$

這會平滑快速變化、限制上升時間，並在跳變附近產生過衝、下衝與振鈴；主要振盪尺度約為 $1/f_{\max}=2/F$。由於 sinc 核在正、負時間都有尾巴，即使原系統因果，重建結果也可能出現到達時間之前的微小響應，這是有限頻寬造成的表觀非因果現象。[L02, pp.10–11]

若截斷前頻譜已非常小，誤差會較輕；否則應增加可靠的量測頻寬或採用適當窗函數。頻寬截斷造成的振鈴，與頻率格距過大導致的時間週期回繞是兩個不同問題：前者主要受 $f_{\max}$ 影響，後者主要受 $\Delta f$ 影響。[L02, pp.5, 10–12]

<a id="q12"></a>

### Q12 What's the purpose of using window function in IDFT?

頻域窗讓頻譜接近截斷邊界時逐漸衰減，減少硬截斷產生的時域振鈴。若 $x_w(t)=\mathcal F^{-1}\{W(f)X(f)\}$，則 $x_w(t)=w(t)*x(t)$；選擇旁瓣衰減較快的 $w(t)$，可以減少遠離主要響應的假振盪。[L02, p.12]

代價是時域主瓣變寬、上升時間變慢，窄脈衝和相近反射可能被混合。通常令 $W(0)=1$，以保持 DC 增益與緩慢變化部分的準確性；對實數波形還要保持 $W(-f)=W^*(f)$，常用實數偶對稱窗。這裡不能不加區分地套用時域頻譜分析的 coherent-gain 補償，因為頻域窗是在有意改變重建波形的高頻內容。[L02, pp.12–15]

<a id="topic-03"></a>

## Sectional convolution

<a id="q13"></a>

### Q13 What's the physical meaning of convolution?

對線性時不變（LTI）系統，$h(t)$ 是單位脈衝輸入造成的輸出。任意輸入可看成許多不同時刻的微小脈衝相加；位於 $\tau$ 的脈衝具有權重 $x(\tau)d\tau$，其輸出便是延遲且縮放的 $x(\tau)h(t-\tau)d\tau$。把所有貢獻相加就得到

$$
y(t)=(x*h)(t)=\int_{-\infty}^{\infty}x(\tau)h(t-\tau)\,d\tau.
$$

這表示現在的輸出通常保留著過去輸入的影響，也就是通道的「記憶」。等價的圖形算法是把一個函數反轉、平移、與另一個相乘後積分；頻域則變成 $Y(f)=X(f)H(f)$。若用離散加總近似連續卷積，通常要乘 $\Delta t$；但如果使用的離散濾波核已包含此因子，就不能再乘一次。[L03, pp.2–4]

<a id="q14"></a>

### Q14 What's the max span of the convolution of two functions?

若 $x(t)$ 只在 $[a,b]$ 非零，$h(t)$ 只在 $[c,d]$ 非零，則卷積的非零範圍包含在 $[a+c,b+d]$。因此兩個支撐長度分別為 $T_1$、$T_2$ 的函數，卷積的最大支撐長度為 $T_1+T_2$；講義的因果例子 $[0,T_1]$ 與 $[0,T_2]$ 對應輸出 $[0,T_1+T_2]$。[L03, p.5]

離散序列則要注意端點計數：長度分別為 $L_x$、$L_h$ 的序列，完整線性卷積長度為 $L_x+L_h-1$。特殊抵消可能讓實際非零範圍更短；若脈衝響應原本有無限長尾巴，任何有限長度都只是截斷近似，必須檢查被捨棄的尾巴是否會影響精度。[L03, pp.5, 13]

<a id="q15"></a>

### Q15 What's sectional convolution and what's the difference to standard convolution?

Sectional convolution 將很長的輸入分成較短區塊，各自與同一個脈衝響應做卷積，再依線性疊加組合。講義使用 overlap-add：每個輸入區塊有 $L$ 點，濾波核有 $P$ 點，則每塊輸出有 $L+P-1$ 點；把它放回原來的時間位置，與相鄰區塊重疊的 $P-1$ 點必須**相加**，不能直接串接。[L03, pp.2, 6]

若每塊採 FFT 計算，轉換長度需滿足 $M\geq L+P-1$，避免 FFT 的循環卷積尾巴繞回區塊前端。若頻率格距已由通道資料固定，則 $M\Delta t=1/\Delta f$ 也必須足以容納區塊與響應尾巴。保持相同取樣、濾波核與邊界處理時，sectional convolution 和一次做完整線性卷積的結果相同，差別是記憶體需求、運算效率與可串流處理的能力。[L03, pp.7, 10–13]

切塊本身不代表每塊都要乘 Hann 等窗，因為這會改變輸入訊號；除非另有能完整重建的加權分段設計。對有限上升時間的 NRZ，區塊第一個轉換也必須接續前一塊的最後電位，不能每塊都重新假設從零開始。通道由有限頻寬資料轉成脈衝響應時所用的頻域窗，則是另一個獨立的模型處理。[L03, pp.6–9; L02, pp.12–15]

<a id="topic-04"></a>

## Spectra of digital signals

<a id="q16"></a>

### Q16 What's the meaning of saying "spectrum of random signal"?

隨機訊號沒有唯一固定的位元序列；每次實現（realization）所得到的波形與有限長 FFT 都會不同。因此談「隨機訊號的頻譜」，必須先指定機率模型，例如各位元獨立、0 與 1 等機率，再討論各種實現的**頻譜功率平均**。不能直接平均複數 FFT 當成 PSD，因為不同實現的相位可能互相抵消，即使平均 FFT 為零，訊號功率仍然存在。[L04, pp.2–4]

持續存在的 NRZ 通常是有限平均功率、無限總能量的訊號，適合用 PSD：$S_x(f)=\lim_{T\to\infty}E\{|X_T(f)|^2\}/T$。固定符元時鐘的 PAM／NRZ 嚴格說可具有週期平穩特性；課堂常用的 PSD 可理解為對時鐘相位做時間平均，或把時間原點在一個 UI 內均勻隨機化後的結果。編碼或重複圖樣會改變統計相關性，不能一律套用獨立隨機位元模型。[L04, pp.2–5, 18–21]

<a id="q17"></a>

### Q17 What's the meaning of power spectral density (PSD)?

PSD 描述平均平方量如何分布在頻率上。若 $x(t)$ 是電壓，本篇的雙邊 $S_x(f)$ 單位是 $\mathrm{V^2/Hz}$，滿足

$$
\int_{-\infty}^{\infty}S_x(f)\,df
=E\{x^2(t)\}=V_{\mathrm{rms}}^2.
$$

若這是實際加在電阻 $R$ 兩端的電壓，功率 PSD 才是 $S_x(f)/R$，單位 $\mathrm{W/Hz}$；講義直接使用 W／Hz 的情況，可視為已採用 $1\ \Omega$ 正規化。含非零平均值時，PSD 可包含 DC 的 delta 線，其積分貢獻平均值平方，不能只把連續部分積分後就稱為總功率。[L04, pp.3–4, 21, 25]

對廣義平穩訊號，Wiener–Khinchin 關係給出 $S_x(f)=\mathcal F\{R_x(\tau)\}$。經過 LTI 系統後，$S_y(f)=|H(f)|^2S_x(f)$，所以只要知道輸入 PSD 和系統傳輸函數，就能算出輸出平均功率分布。PSD 與有限能量脈衝的 energy spectral density $|X(f)|^2$ 不同，後者積分得到總能量而非平均功率。[L04, pp.3–4, 25–26]

<a id="q18"></a>

### Q18 What is the theoretical PSD of random NRZ signal?

先定義 $x(t)=\sum_n a_n p(t-nT_b)$，其中 $p(t)$ 是單位振幅脈衝，$a_n$ 是以伏特表示的獨立同分布電位，平均值為 $\mu_a$、變異數為 $\sigma_a^2$。時間平均的雙邊 PSD 一般為

$$
S_x(f)=\frac{\sigma_a^2}{T_b}|P(f)|^2
+\frac{\mu_a^2}{T_b^2}
\sum_{k=-\infty}^{\infty}|P(k/T_b)|^2
\delta\!\left(f-\frac{k}{T_b}\right).
$$

第一項是連續隨機部分，第二項是平均波形的離散譜線。對 0–100% 線性上升／下降時間同為 $T_r$ 的 NRZ 梯形脈衝，$P(f)=T_b\operatorname{sinc}(fT_b)\operatorname{sinc}(fT_r)e^{-j\pi f(T_b+T_r)}$。若電位是等機率的 $\pm V_{pp}/2$，則 $\mu_a=0$、$\sigma_a^2=V_{pp}^2/4$，因此

$$
\boxed{S_{\mathrm{NRZ}}(f)=\frac{V_{pp}^2T_b}{4}
\operatorname{sinc}^2(fT_b)\operatorname{sinc}^2(fT_r)}.
$$

若電位改為等機率的 $0$ 與 $V_{pp}$，連續部分相同，另加 $(V_{pp}^2/4)\delta(f)$；本脈衝在非零 $k/T_b$ 處的零點使其餘平均波形譜線消失。若 0／1 機率不等，應重新代入平均值與變異數；若位元相關或經過編碼，則需用符元自相關求 PSD。以上是電壓 PSD，求負載功率時再除以電阻。[L01, p.7; L04, pp.4–7, 18, 21]

<a id="q19"></a>

### Q19 How to compute the PSD given timing waveform?

對取樣率 $f_s$ 的電壓波形，取長度 $N_0$ 的資料區塊，加窗 $w[n]$，再求 $M\geq N_0$ 點 FFT：$F_r[k]=\operatorname{FFT}_M\{w[n]x_r[n]\}$。將 $K$ 個區塊或實現的**線性功率**平均，可得雙邊估計式

$$
\widehat S_{x,2}(f_k)
=\frac{1}{Kf_s\sum_{n=0}^{N_0-1}w^2[n]}
\sum_{r=1}^{K}|F_r[k]|^2,
\qquad \Delta f=\frac{f_s}{M}.
$$

這等同於先用 $\Delta t\,F_r[k]$ 近似 FT，再除以記錄時間 $T_0=N_0/f_s$ 與窗的 power gain。分母使用原始窗的平方加總，不能改成補零後的點數，也不能用 coherent gain 的平方取代。先平均線性功率再轉 dB，與直接平均 dB 數值不同；只有當目標是去除平均值後的 AC PSD，才應先扣除各段平均值。[L04, pp.6–7, 15–16]

實數訊號轉為單邊 PSD 時，正頻率部分乘 2，但 DC 與偶數長度時的 Nyquist 點不乘 2。核對正規化時，單段的離散頻率積分應等於窗加權平均平方量 $\sum_n w^2[n]x^2[n]/\sum_n w^2[n]$；矩形窗時才直接變成未加權平均平方量。增加區塊平均次數可降低估計變異；延長每個區塊可改善解析度，兩者在固定總資料量下有所取捨。重疊區塊並非完全獨立，實際變異縮減不會總是恰好 $1/K$；補零則只讓頻率取樣更密。[L04, pp.15–16, 37–43]

<a id="q20"></a>

### Q20 How to correlate the computed PSD with spectrum analyzer measurements?

頻譜分析儀通常顯示某個解析濾波器通帶內的功率，而非直接顯示 $\mathrm{V^2/Hz}$。若雙邊電壓 PSD $S_{v,2}(f)$ 指的是實際負載 $R$ 上的電壓，且在正頻率 $f_c$ 附近夠平坦，則

$$
P_{\mathrm{SA}}(f_c)\approx
\frac{2S_{v,2}(f_c)B_{\mathrm{ENBW}}}{R}
=\frac{S_{v,1}(f_c)B_{\mathrm{ENBW}}}{R}.
$$

因子 2 合併實數訊號的正、負頻率功率，這個近似不直接套用在 DC。更一般的做法是用分析儀實際濾波器的 $|H_{\mathrm{RBW}}(f)|^2$ 對功率 PSD 積分。講義以 RBW 近似頻寬；精確比較雜訊／寬頻功率時應使用 ENBW，因為傳統掃頻分析儀標示的 RBW 可能是 3 dB 頻寬，與 ENBW 不一定相同。[L04, pp.8, 25; W01]

換算單位時，$P_{\mathrm{dBm}}=10\log_{10}(P_{\mathrm{SA}}/1\,\mathrm{mW})$，而 $V_{\mathrm{rms}}=\sqrt{P_{\mathrm{SA}}R}$。若模擬使用的是訊號源開路電壓，必須先計入源阻抗與負載分壓；匹配時負載電壓為開路電壓的一半，功率因而是直接拿開路電壓平方除以 $R$ 所得值的四分之一。對落在通帶內的單一離散譜線，應比較其**積分線功率**，不能把線峰值當作平坦 PSD 再任意乘 RBW；還要對齊儀器的平均方式、檢波器與濾波器設定。[L02, p.8; L04, pp.8, 14]

<a id="q21"></a>

### Q21 What is the bandwidth of NRZ signals?

NRZ 一般沒有「超過某頻率後完全為零」的有限頻寬，必須先說明使用哪種定義。對零上升時間、獨立等機率 NRZ，其連續 PSD 為 sinc 平方形，第一個零點在 $f=1/T_b=R_b$；約 90% 的 AC 功率位於雙邊 $[-R_b,R_b]$，約 95% 位於 $[-2R_b,2R_b]$。對單邊 PSD，等價地說是從 0 積分到 $R_b$ 或 $2R_b$。[L04, pp.5, 9]

有限上升時間再乘上 $\operatorname{sinc}^2(fT_r)$，使高頻衰減更快；工程上也常用 $1/(\pi T_r)\approx0.3/T_r$ 表示高頻包絡開始更快下降的轉折尺度，但它與第一零點、功率涵蓋頻寬並不相同，且必須確認 $T_r$ 的定義。講義建議通道資料至少量到約 $R_b$、最好到 $2R_b$，這是做 NRZ 波形模擬的頻寬涵蓋準則，不等於所有通道都必須具有同樣大的 3 dB 頻寬。[L04, pp.9, 12]

<a id="q22"></a>

### Q22 What's the difference between the spectrum of random NRZ and clock?

獨立隨機 NRZ 的功率主要分布在**連續 PSD**；理想時鐘是週期訊號，功率集中在基頻及諧波的**離散譜線**。若把 `1010...` 視為位元間隔 $T_b$ 的 NRZ，完整週期是 $2T_b$，所以時鐘基頻為 $f_0=1/(2T_b)=R_b/2$，不是 $R_b$。50% duty cycle、對稱邊緣的時鐘除了 DC 外只有奇次諧波，偶次諧波為零。[L04, pp.2, 11]

若 $x(t)=\sum_n c_ne^{j2\pi nf_0t}$，其雙邊 PSD 為 $\sum_n|c_n|^2\delta(f-nf_0)$；對實數時鐘，第 $n>0$ 諧波在電阻 $R$ 上的單邊線功率是 $2|c_n|^2/R$。相同峰對峰振幅與線性邊緣下，諧波包絡同樣受到 $\operatorname{sinc}(fT_b)\operatorname{sinc}(fT_r)$ 控制，但功率集中於窄譜線，所以常對 EMI 峰值特別重要；寬頻隨機訊號在每個 RBW 內的功率則約隨 RBW 增加。[L04, pp.11–14]

有限長 PRBS 若被週期重複，也會產生間距為「圖樣重複頻率」的離散譜線。當圖樣很長、線距比儀器 RBW 小時，量測外觀可能接近連續隨機 PSD；因此比較時要交代圖樣長度、是否重複及量測解析頻寬。[L04, pp.2, 7–8, 14]

<a id="topic-05"></a>

## Eye diagram basics

<a id="q23"></a>

### Q23 How is eye diagram constructed in simulation?

先產生足夠長的資料圖樣，經過通道與等化器得到接收波形，排除啟動暫態，再以已知 UI $T_b$ 和參考時刻 $t_{\mathrm{ref}}$ 切成 1 UI 或 2 UI 片段，將所有片段的相對時間對齊後疊畫。固定速率眼圖可用 $\tau=(t-t_{\mathrm{ref}})\bmod T_b$ 摺疊時間；若要看 2 UI，可用相鄰視窗或在統計穩定的條件下重複一個 UI 的統計結果。[L05, pp.4, 14, 17]

統計眼圖則把相對時間與電壓分箱，累計每個格子的通過次數。每個時間欄除以該欄的總樣本數，可得到「每個電壓 bin 的機率」；若要稱為以 $\mathrm{V^{-1}}$ 為單位的電壓機率密度，還需除以 bin 寬 $\Delta V$。圖樣長度、時間步距、電壓分箱與觸發相位都會影響圖像，稀少的錯誤不能因顏色太淡或分箱太粗就忽略。[L05, pp.17–18]

<a id="q24"></a>

### Q24 How is eye diagram constructed in reality?

示波器以與接收端相關的時鐘作為時間參考：可以使用額外傳送的 forwarded clock，也可以從資料的轉換邊緣以 clock-data recovery（CDR）恢復時鐘。以時鐘通過指定門檻的時刻觸發，擷取資料片段後疊畫，就形成 triggered eye；實作可由即時擷取的記錄處理，或由多次觸發取樣累積。[L05, pp.15–17]

若發送端頻率有微小誤差、漂移或 SSC，硬用一個不準確的固定 UI 長時間摺疊，會造成相位不斷累積，看起來眼圖關閉，即使真正接收器能跟隨時鐘而正常工作。量測需交代觸發／CDR 的追蹤特性、參考位置及儀器頻寬，因為眼圖表示的是**資料相對於所選時鐘**的品質，不是與時間參考無關的唯一圖形。[L05, pp.14–16]

<a id="q25"></a>

### Q25 Basic eye metric: EH, EW, jitter

**Eye height（EH）**是選定取樣時間上的垂直開口；以有限記錄的最小開口定義，可寫成 $EH=\min(v\mid1)-\max(v\mid0)$，兩個集合都取自同一相位。它代表可容忍的電壓擾動空間；若結果不大於零，表示兩個電位集合在該相位已重疊。必須分清楚最小 EH、3σ EH 與某個 BER 等高線上的 EH。[L05, pp.5, 18–20]

**Eye width（EW）**是在指定判決電壓上的水平開口，通常是左側轉換最晚侵入的位置與右側轉換最早侵入的位置之間的距離，代表可用的取樣時間範圍。若要求一個具有指定高度的矩形能放進眼圖，則是 box eye-width，不能只量水平中心線上的 EW。[L05, pp.5, 7]

**Jitter** 是轉換或事件時間相對於參考時間的偏移，可用秒、ps 或 UI 表示；峰對峰 jitter 是該偏移的範圍，扣除平均時間偏移後的 RMS jitter 則是標準差。對具有無界高斯尾巴的隨機 jitter，有限記錄的峰對峰值會隨觀察量增加，不能當成絕對上限。只有在相同交越定義與適當對稱假設下，才能使用 $EW\approx UI-J_{pp}$ 這類簡化關係。[L05, pp.5, 18–20]

<a id="q26"></a>

### Q26 What is the 3-sigma definition of EH and EW?

若取樣相位的邏輯 1、0 電壓分布近似高斯，平均值分別為 $\mu_1,\mu_0$，標準差為 $\sigma_1,\sigma_0$，則

$$
EH_{3\sigma}=(\mu_1-3\sigma_1)-(\mu_0+3\sigma_0).
$$

若左右相鄰交越時間的平均值是 $\mu_L,\mu_R$，標準差為 $\sigma_L,\sigma_R$，則 $EW_{3\sigma}=(\mu_R-3\sigma_R)-(\mu_L+3\sigma_L)$。當兩側標準差同為 $\sigma$，且平均間距為一個 UI，才可簡化成 $EW_{3\sigma}=UI-6\sigma$。[L05, p.19]

這種定義使用分布參數，通常比有限記錄的極端值穩定；但「3σ」不是任何波形都通用的安全範圍。對單一高斯分布，$\pm3\sigma$ 約涵蓋 99.73%，單尾約為 0.135%，因此不能把 3σ 開口直接宣稱成 $10^{-12}$ 等級的 BER 保證。若分布多峰或明顯非高斯，應使用合適的混合模型、實際分位數或指定 BER 準則。[L05, pp.19–20]

<a id="topic-06"></a>

## Superposition, ISI, and PDA

<a id="q27"></a>

### Q27 Why is the single-bit response (SBR) important?

對固定 UI、固定發射脈衝形狀與 LTI 通道，單一位元脈衝的輸出 $s(t)$ 就是 single-bit response。此處將 $a_n$ 定為無因次符元（例如 0/1 或 $\pm1$），$p(t)$ 已包含參考電壓振幅。若輸入寫成 $x(t)=\sum_n a_np(t-nT_b)$，則

$$
s(t)=p(t)*h(t),\qquad
y(t)=\sum_n a_ns(t-nT_b).
$$

因此只要模擬一次短的 SBR，就能用平移、縮放與相加預測任意長位元序列的接收波形，並進一步分析 ISI、設計等化器與預測最差眼圖。[L06, pp.2–3, 7]

SBR 與脈衝響應 $h(t)$ 不同：它已經包含發射位元的寬度、振幅與有限上升時間。因此說 SBR「完整描述通道」，應理解為完整描述**所選脈衝族與線性模型下的資料傳輸**；改變 bit rate 或發射脈衝形狀時，通常需要重新計算 SBR。[L06, pp.2–3; L01, pp.5–7]

<a id="q28"></a>

### Q28 What information can be deduced from the SBR?

從 SBR 的到達時間與主要峰值，可估計通道延遲、主要訊號振幅與合理的取樣相位；每隔一個 UI 取樣，則得到 main cursor、pre-cursors 與 post-cursors。主游標以外的振幅與拖尾長度，顯示其他位元能對目前判決造成多少影響；振鈴與延遲尾巴也可反映通道的頻率選擇性和反射效應，但僅憑一個波形不能唯一判定實體缺陷位置。[L06, pp.5–8]

利用這些游標可重建任意位元序列、推算最差眼圖邊界、評估 ISI，並作為 FFE／DFE 係數設計的基礎。最大 SBR 峰值常是取樣相位的起點，但最大峰值不一定同時給出最大的眼高、眼寬或等化後裕量，還需考慮其他游標與判決條件。[L06, pp.7–8, 12–17]

<a id="q29"></a>

### Q29 What is intersymbol inteference (ISI)?

ISI（intersymbol interference）是某個符元的通道響應延伸到其他符元的判決時間，讓接收電壓除了取決於目前資料，也取決於周圍資料。令 $c_i=s(t_s+iT_b)$，在第 $n$ 個符元的取樣時間可寫成

$$
y_n=a_nc_0+\sum_{i\neq0}a_{n-i}c_i.
$$

第一項是目前符元的 main cursor，後面的總和就是 ISI。$i>0$ 表示過去符元留下的 post-cursor ISI；$i<0$ 是相對於主游標的 pre-cursor ISI。這裡的前／後是以主游標及判決時間為基準，不表示物理通道違反因果性。[L06, pp.4–6, 8]

通道損耗使脈衝展寬，反射與頻率選擇性造成拖尾或振鈴，都可能增加 ISI。結果是同樣的邏輯 1 或 0 會因鄰近圖樣不同而落在不同電壓，造成眼圖閉合與資料相依的轉換時間偏移；資料率越高、UI 越短，固定長度的通道記憶通常會跨越更多符元。[L06, pp.4–6, 19–22]

<a id="q30"></a>

### Q30 What is the peak distortion analysis (PDA) and why is it useful?

PDA 用 SBR 的線性疊加，直接找出每個相位最不利的鄰近位元組合，得到 worst-case-1 的最低電壓與 worst-case-0 的最高電壓。對 $a_n\in\{0,1\}$、SBR 已包含邏輯 1 脈衝振幅的模型，令 $c_0(\tau)$ 為目前位元的貢獻、$c_i(\tau)$ 為其他位元的貢獻，則

$$
WC_1(\tau)=c_0(\tau)+\sum_{i\neq0}\min\{0,c_i(\tau)\},
\qquad
WC_0(\tau)=\sum_{i\neq0}\max\{0,c_i(\tau)\}.
$$

因此該相位的最差垂直開口為 $c_0(\tau)-\sum_{i\neq0}|c_i(\tau)|$。直覺上，要壓低目前的 1，就把會造成負貢獻的鄰近位元設為 1；要抬高目前的 0，就選會造成正貢獻的鄰近位元。若改採 $a_n\in\{-1,+1\}$，則兩條邊界為 $c_0-\sum_{i\neq0}|c_i|$ 與 $-c_0+\sum_{i\neq0}|c_i|$，開口變成 $2(c_0-\sum_{i\neq0}|c_i|)$；兩種編碼的係數不能混用。[L06, pp.7–11]

沿整個 UI 重複上述計算，就能得到眼圖內緣、最差 EH，並依判決門檻或指定矩形高度求 EW。它只需一條足夠長的 SBR，不必等待很長的 PRBS 恰好出現罕見的最差圖樣，因此非常適合快速比較線性通道。各相位的最差圖樣可能不同，所以整條 PDA 邊界是一個逐點包絡，不一定由同一段實際資料波形形成。[L06, pp.12–17]

<a id="q31"></a>

### Q31 What are the limitations of PDA?

基本 PDA 給出的是**所建模型內、資料圖樣造成的最差 ISI 包絡**，沒有提供圖樣出現機率，因此不能直接得到統計眼圖、3σ EH／EW 或 BER。基本形式也未包含隨機電壓雜訊、隨機／週期 jitter 與獨立侵擾源的串擾；要納入這些效應，需額外的模型及擴充算法。[L06, p.17]

它依賴固定脈衝的線性疊加，因而假設通道與所包含的收發器模型為線性、時不變，且上升與下降形狀能由同一個脈衝模型描述。若有上升／下降不對稱、飽和、資料相依的驅動器行為，或 DFE 的錯誤傳播等非線性記憶，單一 SBR 的基本 PDA 不足以描述實際系統。[L06, pp.17–18]

此外，逐點選出的最差圖樣若違反實際編碼限制，可能過度保守；反過來，若 SBR 截得太短、忽略長尾或使用不正確的通道模型，則可能把最差開口估得過大。有限 PRBS 模擬的眼圖比 PDA 大，並不自動代表 PDA 算錯，可能只是該記錄尚未遇到允許的最差圖樣。[L06, pp.7–8, 17; L04, p.18]

<a id="topic-07"></a>

## Equalizers

<a id="q32"></a>

### Q32 How to estimate the risetime of a channel based on its S21?

先把通道看成低通系統，找相對於低頻增益下降 $3\,\mathrm{dB}$ 的頻率 $f_{3\mathrm{dB}}$。若能用單極點 $H(s)=H(0)/(1+s\tau)$ 近似，理想階躍的輸出為 $H(0)(1-e^{-t/\tau})$，所以

$$
t_{10\text{--}90}=\tau\ln9\approx2.197\tau,
\qquad f_{3\mathrm{dB}}=\frac1{2\pi\tau},
\qquad t_{10\text{--}90}\approx\frac{0.35}{f_{3\mathrm{dB}}}.
$$

這是**單極點或相近低通形狀的估計**，不是所有通道的定律。多極點、反射凹谷、頻率相依群延遲及振鈴都可能使估計失準；講義例子由頻寬估得約 $95\,\mathrm{ps}$，直接模擬約為 $112\,\mathrm{ps}$。較可靠的做法是保留複數 $S_{21}$ 的振幅與相位，先求脈衝響應、再積分成階躍響應，依最終階躍幅度量測 10% 與 90% 交越。若各埠採共同實數參考阻抗 $Z_0$，且 Thévenin 源串聯電阻與負載均為 $Z_0$，從**理想源的開路電壓**到負載電壓的傳輸函數是 $S_{21}/2$；若輸入變數是已進入通道的入射電壓波，則使用 $S_{21}$。常數尺度不改變正規化上升時間。來源的有限上升時間也要納入，不能把量到的整體上升時間全歸因於通道。[L07, p.2–4]

<a id="q33"></a>

### Q33 Why lossy channel cannot transmit high data rate signals?

資料率增加時，$T_b$ 變短，通道卻仍以自己的頻寬、色散與記憶長度回應。導體損耗、介質損耗通常使高頻比低頻衰減更多，邊緣變慢、單位脈衝拖尾；當相鄰位元的尾巴在取樣時刻重疊，就形成 ISI。由階躍響應 $r(t)$ 可寫出理想矩形位元的 SBR：$p_{\mathrm{out}}(t)=r(t)-r(t-T_b)$。若 $T_b$ 遠小於通道的上升時間，兩個階躍響應高度重疊，主游標變小，其他位元造成的干擾相對變大，眼圖因而閉合。[L07, p.3–5]

因此「不能傳高速」應理解成：在指定振幅、雜訊、等化器與 BER 要求下，存在可用資料率的限制。通道有損耗並不代表一定不能高速傳輸；足夠的接收靈敏度、編碼與等化可以改善結果。反之，等化也不能無限制地補回深凹谷或已低於雜訊的訊號。NRZ 的 Nyquist 頻率為 $R_b/2$，但實際邊緣包含更高頻成分，不能只檢查單一頻率的 $S_{21}$ 就宣判通道可用。[L07, p.5、13–14]

<a id="q34"></a>

### Q34 What is the concept of equalization?

等化是在 Tx、Rx 或兩端加入補償，使**取樣時的主游標相對於 ISI 更突出**，並改善可用眼圖。對線性等化器，理想目標常寫成

$$
H_{\mathrm{eq}}(f)H_{\mathrm{ch}}(f)
\approx K e^{-j2\pi f t_d}
$$

也就是在有用頻帶內，總增益較平坦、相位接近線性，允許一個整體延遲 $t_d$。CTLE 用連續時間的極點、零點補償；FFE 用多個延遲且加權的副本降低游標 ISI；DFE 則根據已判決的過去位元，取消可預測的後游標。DFE 含 slicer，整個系統不是一般 LTI 濾波器，不能為它直接指定一個與輸入無關的 $H_{\mathrm{eq}}(f)$。[L07, p.5–6；L08, p.2–3；L09, p.4–6]

工程目標也不是不計代價地實現 $1/H_{\mathrm{ch}}$。若補償會放大雜訊、超過 Tx 電壓限制、需要太多 taps，或造成更多其他 ISI，就必須折衷；最終應以指定時鐘下的 EH、EW 與錯誤表現驗證。[L07, p.13–14；L08, p.22；L09, p.13–14]

<a id="q35"></a>

### Q35 What is CTLE and how is it realized in circuit?

CTLE 是 continuous-time linear equalizer，以連續時間類比濾波器，補償通道高頻相對衰減較多的特性。典型的一零點、兩極點形式可寫為

$$
H_{\mathrm{CTLE}}(s)
=A_{\mathrm{DC}}
\frac{1+s/\omega_z}{(1+s/\omega_{p1})(1+s/\omega_{p2})},
\qquad \omega_z<\omega_{p1}<\omega_{p2}.
$$

零點與第一極點之間形成相對高頻提升，第二極點限制極高頻的增益。被動 CTLE 可由 RC、RL、RLC 網路實現，主要靠降低低頻，使訊號頻寬相對於 DC 增大；它不需要主動供電來提供功率增益。主動 CTLE 常用帶電阻、電容源極退化的差動放大器：低頻退化較強，高頻時電容降低退化阻抗，因此得到 peaking。實際電路仍受負載電容、功耗、雜訊與線性範圍限制。[L07, p.6–9、11]

設計時應檢查**真正的最高增益**，不能把傳輸式中的某個分子係數直接當成峰值。CTLE 也會處理進入 Rx 的雜訊與串擾，過度提升可能使它們更嚴重；少數極點、零點通常也無法精確補平由多重反射形成的複雜凹谷。[L07, p.8、13–14；L09, p.2–3]

<a id="q36"></a>

### Q36 What is FFE and how is it realized in circuit?

FFE 是 feedforward equalizer，本質為 FIR：把資料的多個延遲副本加權相加。若 tap 間隔 $T$ 通常取為 $1\,\mathrm{UI}$，以主游標標成 $n=0$，則

$$
y_k=\sum_n c_n x_{k-n},
\qquad y(t)=\sum_n c_n x(t-nT).
$$

$n<0$ 是 precursor taps，$n>0$ 是 postcursor taps。這些名稱是**相對主游標**的時間位置；實體電路可加入共同延遲後，使所有分支皆因果，不需要真的預知尚未收到的資料。[L08, p.3、17–18]

Tx 常以移位暫存器產生相鄰位元，再讓多組具有正、負權重的電流型或分段驅動器，在輸出端相加；tap 權重可由驅動電流或開啟的單元數實現。Rx 也可做 FFE，例如以取樣後的多 tap 濾波器實現，但高速類比延遲、雜訊放大與適應性調整的代價不同。FFE 的作用是改變送入通道或接收端的**類比波形**，不只是事後更改數位資料。[L08, p.15、22]

<a id="q37"></a>

### Q37 What are the meaning of pre-emphasis and de-emphasis?

兩者都使轉換附近的高頻成分，相對於長串相同位元的低頻成分更大。以轉換時最大擺幅 $V_{\mathrm{swing}}$ 與長串資料的穩態電平差 $V_{\mathrm{shelf}}$ 表示：pre-emphasis 保留原本穩態電平、增加轉換擺幅；de-emphasis 保留允許的最大擺幅、降低穩態電平。它們可能有相同的相對等化形狀，但絕對電壓尺度不同；電源與可靠度限制使 Tx 常採 de-emphasis。[L08, p.4–6、10]

例如輸入為 $0/1\,\mathrm V$，兩 tap $y_k=c_0x_k+c_1x_{k-1}$。若 $(c_0,c_1)=(1.2,-0.2)$，長串 1 為 $1.0\,\mathrm V$，最大、最小值為 $1.2$、$-0.2\,\mathrm V$，總擺幅為 $1.4\,\mathrm V$；若改成 $(0.8,-0.2)$，長串 1 為 $0.6\,\mathrm V$，最大、最小值為 $0.8$、$-0.2\,\mathrm V$，總擺幅維持 $1.0\,\mathrm V$。後者的 de-emphasis 比率是 $V_{\mathrm{shelf}}/V_{\mathrm{swing}}=0.6$，若用幅度 dB 表示則為 $20\log_{10}0.6\approx-4.44\,\mathrm{dB}$；有些文件以正值稱呼「衰減量」，須先看定義。[L08, p.4–9]

<a id="q38"></a>

### Q38 What is the equivalent transfer function of an FFE?

由延遲的 Fourier transform 性質，FFE 的等效頻率響應為

$$
H_{\mathrm{FFE}}(f)=\sum_n c_n e^{-j2\pi fnT},
\qquad H_{\mathrm{FFE}}(z)=\sum_n c_n z^{-n}.
$$

因此等化後 SBR 是各延遲 SBR 的加權和；對 LTI 通道，Tx FFE 與通道的串接可在頻域相乘。若為了實現 precursor taps 而加入共同延遲，只會另外乘上 $e^{-j2\pi f t_d}$，不改變增益大小。理想 tap-delay 模型的頻率響應以 $1/T$ 為週期，實際驅動器的有限頻寬與邊緣則還需另外包含。[L08, p.3、11–12、17]

對 $(c_0,c_1)=(0.8,-0.2)$、$T=T_b$，$H(0)=0.6$，而在 NRZ Nyquist 頻率 $f=1/(2T_b)$，$H=0.8-(-0.2)=1$，呈現相對高頻提升。注意 $H(1/T_b)=H(0)=0.6$；「位元率」與「Nyquist 頻率」不能混用。[L08, p.11]

<a id="q39"></a>

### Q39 What is the maximum-swing constraint of FFE?

對任意可獨立取 $0$ 或 $A$ 的輸入位元，$y_k=\sum_n c_nx_{k-n}$。最大輸出由正 taps 選 $A$、負 taps 選 0 得到；最小輸出反過來，因此

$$
y_{\max}=A\sum_{c_n>0}c_n,\qquad
y_{\min}=A\sum_{c_n<0}c_n,\qquad
V_{\mathrm{swing}}=A\sum_n|c_n|.
$$

若希望 FFE 後仍不超過原來的峰對峰電壓 $A$，就要求 $\sum_n|c_n|\le1$；充分使用這個上限時，通常正規化成 $\sum_n|c_n|=1$。這是 **L1 norm**，不是帶正負號的 $\sum_n c_n=1$。後者決定 DC gain；在常見負旁 tap 的 de-emphasis 設計中通常小於 1。[L08, p.7–9、17、21]

若輸入是 bipolar $\pm A$，最大、最小值為 $\pm A\sum_n|c_n|$，所以峰對峰電壓為 $2A\sum_n|c_n|$，同樣得到 L1 限制。差動訊號還要分清「每條單端 driver 的 Vpp」與「兩線相減後的差動 Vpp」：互補驅動時後者是前者的兩倍。有限測試資料未碰到最壞組合，不能取代這個任意位元組合的上限。[L10, p.20–21、26–28]

<a id="q40"></a>

### Q40 How to design the FFE tap coefficients?

先取得欲補償系統的 SBR，例如 channel 加 CTLE，並在選定的取樣相位取出 $a_j=p(t_0+jT_b)$。經 FFE 後，游標為離散摺積

$$
r_j=\sum_n c_n a_{j-n}.
$$

Zero-forcing 選擇有限個 $j$，令主游標等於指定值 $D$、其他目標游標為零，得到 $\mathbf A\mathbf c=D\mathbf e_0$，其中 $A_{j,n}=a_{j-n}$。應求解這個線性方程，而非把「把每個 postcursor 各除一次」當作多 tap 的通解。解出後按最大擺幅正規化 $\mathbf c\leftarrow\mathbf c/\sum|c_n|$；正規化保留目標零點，但會降低主游標的絕對幅度。[L08, p.12、17–21]

有限 taps 只能限制部分游標，ZF 可能犧牲太多 Tx 振幅、把其他游標變大，或對通道估測誤差敏感。也可做加權最小平方、含雜訊的 MMSE、訓練序列適應調整，或直接在 tap／電壓限制下最佳化 EH、EW。若後面有 DFE，FFE 往往優先處理 precursor 與重要的短期 ISI，讓 DFE 接手 postcursor，通常比強迫 FFE 單獨歸零所有尾巴更合理。[L08, p.14、22、24；L09, p.14、16]

<a id="q41"></a>

### Q41 What is DFE and how is it realized in circuit?

DFE 是 decision feedback equalizer。它用過去的**接收端判決**估計先前位元對現在造成多少 postcursor ISI，從接收波形扣除後，再判斷目前的位元。對 bipolar 資料，可寫成

$$
z_k=x_k-\sum_{j=1}^{N_d}d_j\hat b_{k-j},
\qquad
\hat b_k=\begin{cases}+1,&z_k\ge V_t,\\-1,&z_k<V_t.\end{cases}
$$

硬體包含 slicer／比較器、用 flip-flops 儲存與延遲判決的回授路徑、tap 權重 DAC 或加權電流單元，以及類比加總器。slicer 使回授值回到離散邏輯電平，不把每次接收的類比雜訊原樣再循環；代價是回授必須在很短的 UI 內完成，時鐘、延遲與建立時間都很關鍵。[L09, p.4–6、13]

「DFE 眼圖」應看 slicer **之前的加總器輸出** $z(t)$，不能只畫已被強制變成理想 0/1 或 $\pm1$ 的判決輸出。實際 DFE 是具記憶的非線性系統，也不能把已知 Tx bits 直接餵入回授，再當成真實 DFE 的結果。[L09, p.4–6、9–10]

<a id="q42"></a>

### Q42 Why do we need DFE in addition to CTLE and FFE?

CTLE 與 FFE 都是線性補償：Rx 提升被衰減的頻帶時，也會提升該處的雜訊與串擾；Tx FFE 雖不直接放大後來加入的 Rx 雜訊，卻受最大擺幅限制，常以降低穩態訊號來換等化。DFE 利用已判決的符號去預測 postcursor，不必用大幅高頻增益來抵消同樣的尾巴，因此可在不連同放大類比輸入雜訊的情況下，降低可預測 ISI。這不等於它能消除接收雜訊，或保證改善任何定義、任何條件下的 SNR。[L08, p.22；L09, p.2–4、13]

典型分工是 CTLE 補償平滑的頻率相依損耗，FFE 處理 precursor 與部分短期 ISI，DFE 清除剩餘 postcursor。DFE 無法用尚未判決的未來位元取消 precursor；若初始眼圖太閉合，錯誤判決會產生錯誤回授，並影響後面位元，形成 error propagation。所以它通常需要前級先提供足以可靠判決的開口，不能單獨拯救所有嚴重失真的通道。[L09, p.9–14]

<a id="q43"></a>

### Q43 How to design the DFE tap coefficients?

固定前級 CTLE、FFE 與接收取樣時鐘後，從**整個 pre-DFE 系統**的 SBR 取得 $a_j=p(t_0+jT_b)$。若資料符號正規化為 $b_k\in\{-1,+1\}$，且

$$
x_k=a_0b_k+\sum_{j\ge1}a_jb_{k-j}+\text{precursor ISI}+n_k,
\qquad z_k=x_k-\sum_{j=1}^{N_d}d_j\hat b_{k-j},
$$

則在過去判決正確時選 $d_j=a_j$，即可取消所涵蓋的 postcursor。若加總器寫成 $z_k=x_k+\sum c_j\hat b_{k-j}$，相同設計應填 $c_j=-a_j$。只有當整個 $x_k$ 先除以 $a_0$ 正規化，才使用無因次 $d_j=a_j/a_0$；直接在伏特域扣除時，taps 的單位也是伏特。[L09, p.5–7；L10, p.30–32]

對 0/1 slicer，也應以「單位 1 符號」對應的 SBR 取得係數，且判為 0 時該次回授為零；對 $\pm1$ slicer，判為負符號仍須提供反方向回授。若量到的是 bipolar 兩電平差所產生的增量脈衝，須先除以 2，才是 $\pm1$ 單位符號的係數；這個因子與匹配電壓源的 $S_{21}/2$ 分壓是兩件不同的事。最後用**遞迴判決**驗證錯誤傳播與眼圖，並檢查回授保持區間及時鐘相位，不能只把 SBR 尾端改成零就宣稱 DFE 已驗證。[L09, p.5–6、9–10；L10, p.25、31]

<a id="topic-08"></a>

## Differential signaling (1)

<a id="q44"></a>

### Q44 What are the pros and cons of differential signaling vs single-end?

差動用兩條線傳一個符號，接收端判讀兩線之差。好處是可抑制兩線上相近的共模干擾、降低對 Tx/Rx 共模 DC 偏壓差的敏感度，並在走線靠近且平衡時，減少對外輻射與部分回流不連續的影響。以兩線各有相同訊號振幅、相同且互不相關的雜訊作比較，訊號相減後加倍，雜訊功率則加倍，所以相對單條線的**功率 SNR 是兩倍，即約 $3\,\mathrm{dB}$**；若雜訊主要為共模，改善可能更大，若比較條件改為相同總功耗或相同總擺幅，就不能直接沿用這個數字。[L10, p.7、9]

代價是兩條走線、兩套驅動與接收相關電路、更多空間及較嚴格的配對要求。長度、負載、過孔或材料不對稱會產生 skew 與模式轉換，削弱共模抑制並增加 EMI；差動也仍然需要良好的回流路徑。AC 耦合可隔開 DC 偏壓，但差動本身不保證消除長連續碼造成的 droop／baseline wander，仍需考慮編碼、耦合電容與資料圖樣。[L10, p.2–9、12–13]

<a id="q45"></a>

### Q45 What are differential-mode (DM) and common-mode (CM)?

對同一截面的兩條線，DM 是大小相等、極性相反的分量；CM 是大小相等、極性相同的分量。任意兩線電壓皆可唯一分解成

$$
\begin{bmatrix}V_1\\V_2\end{bmatrix}
=V_c\begin{bmatrix}1\\1\end{bmatrix}
+\frac{V_d}{2}\begin{bmatrix}1\\-1\end{bmatrix}.
$$

純 DM 的共模電壓為零，純 CM 的差動電壓為零；若還有 DC bias，通常是把 DM 疊加在一個共同偏壓上。理想差動接收器只對 $V_d$ 有反應，但真實接收器的 CMRR、共模輸入範圍及通道對稱性有限。通道不對稱時，DM 可轉成 CM，CM 也可轉成 DM；所以不能只看源端是否送了互補訊號。[L10, p.13–17、21–24]

<a id="q46"></a>

### Q46 What are the definitions of DM/CM voltages/currents?

兩條線的電流採相同縱向參考方向；作為網路埠時，兩者皆採流入網路。功率一致的定義是

$$
V_d=V_1-V_2,\qquad I_d=\frac{I_1-I_2}{2},
\qquad V_c=\frac{V_1+V_2}{2},\qquad I_c=I_1+I_2.
$$

反變換為 $V_{1,2}=V_c\pm V_d/2$、$I_{1,2}=I_c/2\pm I_d$，並滿足 $V_1I_1^*+V_2I_2^*=V_dI_d^*+V_cI_c^*$。係數不能任意刪掉，否則功率與阻抗定義就不一致。對參考平面上兩條各為 $Z_0$ 的不耦合參考線，純 DM 與 CM 的參考阻抗分別為

$$
Z_d=2Z_0,\qquad Z_c=\frac{Z_0}{2}.
$$

因此 $Z_0=50\,\Omega$ 對應 $Z_d=100\,\Omega$、$Z_c=25\,\Omega$。這是 mixed-mode S 參數的參考阻抗；實際耦合線的模態特性阻抗則另由其幾何與材料決定，例如採相同定義時 $Z_{\mathrm{diff}}=2Z_{\mathrm{odd}}$。[L10, p.14–15]

<a id="q47"></a>

### Q47 What is mixed-mode S matrix, and how is it related to single-ended S matrix?

Mixed-mode S matrix 以差動、共模的**功率波**為埠變數，描述同模式的反射、傳輸與模式轉換。本課四埠配對固定為輸入端 $(1,3)$、輸出端 $(2,4)$，模態順序為 $[c_1,d_1,c_2,d_2]$。對各單端埠皆為相同實數 $Z_0$，以及 Q46 的 $Z_d$、$Z_c$，令

$$
\mathbf T=\frac1{\sqrt2}
\begin{bmatrix}
1&0&1&0\\1&0&-1&0\\0&1&0&1\\0&1&0&-1
\end{bmatrix},
\quad
\mathbf a_m=\mathbf T\mathbf a,\quad
\mathbf b_m=\mathbf T\mathbf b.
$$

由 $\mathbf b=\mathbf S\mathbf a$、$\mathbf T^{-1}=\mathbf T^T$ 得 $\boxed{\mathbf S_m=\mathbf T\mathbf S\mathbf T^T}$。$1/\sqrt2$ 是功率波正規化，與電壓定義中的 $1/2$ 不能混為一談。$S_{dd21}$ 表示 DM 輸入到 DM 輸出，$S_{cd21}$ 是 DM 輸入轉 CM 輸出；前一個模式下標是輸出，後一個是輸入。埠排序、配對或參考阻抗改變時，必須跟著改轉換，不可機械套同一個矩陣。[L10, p.15–17]

<a id="q48"></a>

### Q48 How to quickly compute Sdd21, Sdd11, etc from single-ended S parameters?

對 Q47 的埠配置，可用「DM 取差、CM 取和，再乘 $1/2$」快速算出：

$$
\begin{aligned}
S_{dd21}&=\tfrac12(S_{21}-S_{23}-S_{41}+S_{43}),\\
S_{dd11}&=\tfrac12(S_{11}-S_{13}-S_{31}+S_{33}),\\
S_{cc21}&=\tfrac12(S_{21}+S_{23}+S_{41}+S_{43}),\\
S_{dc21}&=\tfrac12(S_{21}+S_{23}-S_{41}-S_{43}),\\
S_{cd21}&=\tfrac12(S_{21}-S_{23}+S_{41}-S_{43}).
\end{aligned}
$$

例如 $S_{dd21}$ 的輸入端用 $(+1,-1)$ 激發埠 1、3，輸出端再取埠 2 減埠 4；兩次功率波的 $1/\sqrt2$ 合起來成為 $1/2$。要對**複數值**相加減，不能先取 dB 或只取絕對值。這些是參考阻抗與功率波已正規化的 S 參數；若由理想 Thévenin 差動電壓源、經匹配串聯電阻驅動匹配通道，源電壓到負載差動電壓還有另一個 $1/2$，即 $H_{\mathrm{source}\to\mathrm{load}}=S_{dd21}/2$。[L10, p.17、20–22]

<a id="q49"></a>

### Q49 What is the SBR of a differential channel?

先把數據寫成 $b_k\in\{-1,+1\}$，源端差動波形為 $v_{d,\mathrm{src}}(t)=\sum_k b_kp(t-kT_b)$，其中 $p(t)$ 是代表「單位正符號」的有限上升時間脈衝。差動 SBR 定義為這個**純 DM 基底脈衝**通過通道後的輸出；在每線源、負載皆匹配時，

$$
p_d(t)=\mathcal F^{-1}\!\left\{\frac{S_{dd21}(f)}2P(f)\right\},
\qquad v_{d,\mathrm{out}}(t)=\sum_k b_kp_d(t-kT_b).
$$

例如 $0<T_r\le T_b$，$p(t)$ 是由 $t=0$ 開始上升、於 $T_r$ 達到 $A$、於 $T_b$ 開始下降且於 $T_b+T_r$ 回到零的梯形，則在本文的 FT、normalized sinc 定義下，$P(f)=AT_b\,\mathrm{sinc}(fT_b)\mathrm{sinc}(fT_r)e^{-j\pi f(T_b+T_r)}$。這裡的「其他符號為零」是線性分解使用的數學基底，並不是聲稱真實 bipolar 資料有第三種資料電平；也不能把一條線送單一 1、另一條永遠送其補數的整段波形，直接當成純單位 DM 脈衝。[L10, p.20–25]

若源含共模輸入且通道有 CM-to-DM 轉換，負載還有 $v_{c,\mathrm{src}}*\mathcal F^{-1}\{S_{dc21}\}$ 的貢獻，須分開處理。當共模偏壓固定且系統已穩定時，資料相關的等化通常以 DM SBR 為主；源端兩路不匹配或共模隨資料變動時，單一 DM SBR 可能不足。[L10, p.22–25、28–30]

<a id="q50"></a>

### Q50 How does the post-FFE waveform look like in a differential channel?

把 bipolar 符號 $b_k=\pm1$ 經 FFE 得 $q_k=\sum_jc_jb_{k-j}$，源端差動波形便是 $v_{d,\mathrm{FFE}}(t)=\sum_kq_kp(t-kT_b)$。它通常有多個振幅電平：轉換附近較大，長串相同位元的 shelf 較小，並可依 precursor taps 出現 preshoot。它仍是類比波形，有限 rise/fall time、過衝與 driver 限制皆應保留。[L10, p.26–28；L08, p.9、17、23]

例如兩 tap $(0.8,-0.2)$，連續相同資料時 $q_k=\pm0.6$；剛發生轉換時 $q_k=\pm1$。若每線 driver 原始 Vpp 為 $A$，源端差動電平可達 $\pm A$，差動 Vpp 為 $2A$；L1 正規化使最壞資料組合仍遵守這個上限。兩條線可寫成 $v_1=v_c+v_d/2$、$v_3=v_c-v_d/2$。如果兩路都對原來的 0/1 波形做同一線性 FFE，共模 DC 還會受到 $\sum c_j$ 縮放；若電路另外控制 bias，則需依實際 bias 電路決定，共模不宜憑空假設。[L10, p.26–28]

<a id="q51"></a>

### Q51 What is the slicer output of a DFE for differential channel?

對本課以 $b_k\in\{-1,+1\}$ 表示的差動資料，DFE slicer 的正規化輸出應為 **$+1$ 或 $-1$**。在偏壓對稱、差動輸出繞零且沒有需補償的 offset 時，常用門檻 $V_t=0$：$z_k\ge0$ 判 $+1$，否則判 $-1$。這個 $\pm1$ 是回授符號的數值，不表示輸入 driver 一定是 $\pm1\,\mathrm V$；實際電壓由 SBR 與 DAC tap 尺度決定。[L10, p.30–32]

因此不能沿用單端 0/1 寫法，把「判為負資料」當成不做任何事。在 $z=x-\sum_jd_j\hat b$ 的約定下，$+1$ 扣掉 $d_j$，$-1$ 則加回 $d_j$。非對稱雜訊、offset 或其他接收設計可能使最佳門檻偏離零；取樣相位與回授延遲也必須和擷取 taps 時一致。[L09, p.5–6；L10, p.31]

<a id="topic-09"></a>

## Crosstalk

<a id="q52"></a>

### Q52 Frequency domain solution of NEXT and FEXT

先固定端口：aggressor 的近、遠端為 1、2，victim 的近、遠端為 3、4，兩線沿同一方向，victim 電壓皆對同一參考地取正。對均勻、互易、對稱耦合線，可用 even／odd 模態分解。若模態特性阻抗為 $Z_e,Z_o$、傳播長度為 $\ell$、各外接埠參考阻抗為 $Z_0$，令

$$
\Gamma_r=\frac{Z_r-Z_0}{Z_r+Z_0},\quad q_r=e^{-\gamma_r\ell},\quad
R_r=\frac{\Gamma_r(1-q_r^2)}{1-\Gamma_r^2q_r^2},\quad
T_r=\frac{(1-\Gamma_r^2)q_r}{1-\Gamma_r^2q_r^2},\quad r=e,o.
$$

則 $\boxed{S_{31}=\tfrac12(R_e-R_o)}$ 是 NEXT，$\boxed{S_{41}=\tfrac12(T_e-T_o)}$ 是 FEXT，另有 $S_{21}=\tfrac12(T_e+T_o)$。這些式子已包含每個模態在兩端失配後的重複反射；不可在未檢查端接條件時，一開始就只留下單次傳播項。[L11, p.6–9]

在弱耦合、小模態反射、低損耗，且 $\Gamma_e\approx-\Gamma_o\approx\Gamma_N$、近似共同單程延遲 $\tau$ 時，NEXT 可近似為 $S_{31}\approx\Gamma_N(1-e^{-j2\omega\tau})$；FEXT 可近似為 $S_{41}\approx\tfrac12(e^{-j\omega\tau_e}-e^{-j\omega\tau_o})$。這裡 $\omega=2\pi f$，$\tau_{e,o}=\ell/v_{e,o}$。S 參數是入射功率波的比例；若問的是匹配理想電壓源到 victim 負載的電壓傳輸，還應使用 $S_{31}/2$ 或 $S_{41}/2$。[L11, p.7、9、18]

<a id="q53"></a>

### Q53 Frequency domain characteristics of NEXT and FEXT

在 Q52 的近似下，

$$
|S_{31}|\approx2|\Gamma_N|\,|\sin(\omega\tau)|.
$$

低頻時 NEXT 約正比於 $f$，在振幅 dB 圖上呈 $20\,\mathrm{dB/decade}$；它不是隨頻率永遠增加，而有相長、相消的波紋。第一個非 DC 零點是 $\omega\tau=\pi$，即 $f=1/(2\tau)$，此時耦合長度約為該模態半波長。損耗、色散、模態延遲不等與端接失配會改變零點深度及位置。[L11, p.7、15–16]

令 $\Delta\tau=\tau_o-\tau_e$、$\bar\tau=(\tau_e+\tau_o)/2$，則

$$
S_{41}\approx j e^{-j\omega\bar\tau}\sin\!\left(\frac{\omega\Delta\tau}{2}\right)
\approx j\omega\frac{\Delta\tau}{2}e^{-j\omega\bar\tau}
\quad(|\omega\Delta\tau|\ll1).
$$

因此 FEXT 在小模態相位差區域也約正比於 $f$，大小取決於**兩模態的延遲差**。若均勻介質中兩模態等速，一階 FEXT 可相消；接頭、不均勻性與較高階失配效應仍可能留下串擾。有限的 $\sin$ 響應也說明微分近似不能外推至任意高頻。注意 $\Delta\tau=\ell(1/v_o-1/v_e)$ 的單位是秒；不能漏掉長度 $\ell$，否則 $j\omega\Delta\tau$ 不再無因次。[L11, p.9、15–17]

<a id="q54"></a>

### Q54 Time domain solution of NEXT and FEXT

令 $x(t)$ 是**已進入 aggressor 近端的入射電壓波形**，不是尚未經匹配分壓的理想源電壓。在 Q52 相同的弱耦合、近似非色散條件下，頻域式反變換得到

$$
v_{\mathrm{NEXT}}(t)\approx\Gamma_N[x(t)-x(t-2\tau)],
\qquad
v_{\mathrm{FEXT}}(t)\approx\frac12[x(t-\tau_e)-x(t-\tau_o)].
$$

再令模態延遲差遠小於輸入變化時間，就有

$$
v_{\mathrm{FEXT}}(t)\approx\frac{\Delta\tau}{2}
\left.\frac{dx(u)}{du}\right|_{u=t-\bar\tau}.
$$

這也直接給出單位檢查：秒乘伏特／秒，結果為伏特。若輸入資料由理想 Thévenin 源經匹配電阻產生，應先令 $x=v_{\mathrm{src}}/2$，不要又在末端把同一分壓重算一次。一般有損、色散或明顯反射時，應回到完整 $H_{\mathrm{NEXT/FEXT}}(f)$，對實際有限邊緣波形做卷積，而非強行套兩個延遲或微分式。[L11, p.7、9、16、18]

<a id="q55"></a>

### Q55 Time domain characteristics of NEXT and FEXT

對一個上升階躍，NEXT 是輸入與延遲 $2\tau$ 的副本之差：耦合區由近端開始時，近端串擾幾乎立即出現；理想零上升時間時形成高度約 $\Gamma_N\Delta V$、持續約 $2\tau$ 的脈衝。有限上升時間會使脈衝邊緣變圓；當 $2\tau\ll T_r$，峰值約隨耦合長度增加；當 $2\tau\gg T_r$，峰值趨近飽和，增加長度主要延長持續時間。若耦合區不在接頭起點，還須加上到達耦合區的傳播延遲。[L11, p.7、16]

FEXT 約在單程延遲 $\bar\tau$ 後到達，主要集中在 aggressor 的轉換附近。微分近似給出線性斜坡的峰值約 $\Delta V\,\Delta\tau/(2T_r)$，因此邊緣越快、模態速度差或耦合長度越大，FEXT 通常越大；上升與下降邊緣的極性相反。若 $|\Delta\tau|$ 已可和 $T_r$ 比較，應用兩個延遲階躍之差，不能繼續用微分近似。資料串流中的 NEXT、FEXT 由許多這種貢獻疊加，且會隨 aggressor 相對 victim 的碼型與相位而改變，不應直接視為獨立高斯雜訊。[L11, p.9、12–17]

<a id="q56"></a>

### Q56 Effect of coupling on insertion loss

即使相鄰線沒有送出 aggressor 資料，耦合仍改變原線的模態阻抗、傳播與能量分配；因此「無 aggressor」不等於「和單獨一條線的 $S_{21}$ 相同」。在小模態反射、近似等幅的例子中，

$$
S_{21}\approx\frac12(e^{-j\omega\tau_e}+e^{-j\omega\tau_o})
=e^{-j\omega\bar\tau}\cos\!\left(\frac{\omega\Delta\tau}{2}\right).
$$

even、odd 到遠端的相位差達 $\pi$ 時，兩者在原線上相消，形成 insertion-loss 凹谷／零點；在弱損耗模型中，同頻率的 FEXT 則相長。較長耦合區或較大的模態延遲差，會使第一個這類凹谷降到較低頻。[L11, p.3–6、8、14]

這種 $|S_{21}|$ 降低不一定全是導體或介質把能量變熱：對理想無損且匹配參考埠的四埠網路，$\sum_{i=1}^4|S_{i1}|^2=1$，能量也可分配到反射及其他線。實際通道常同時包含真正耗散損耗與耦合效應；分析時應看完整多埠資料及端接，而非只替原本單線 $S_{21}$ 額外加一個雜訊源。[L11, p.6、8、14]

<a id="q57"></a>

### Q57 How to correctly simulate the effects of crosstalk?

先建完整多埠通道並保留所有線與端接。對 LTI 系統，可分別計算每個發射源到 victim 的響應再相加：$V_i(f)=\sum_jH_{ij}(f)V_{Tj}(f)$。若每個源與負載都等於實數參考阻抗 $Z_0$，則 $H_{ij}=S_{ij}/2$，包括 victim 自身傳輸與所有 NEXT／FEXT 路徑；未啟動的源應設其理想源電壓為零，仍保留串聯源阻抗，不能任意把該埠開路。若改為差動，須用一致埠配對的 mixed-mode 通道及正確源尺度。[L11, p.18]

失配時不能直接用 $S_{ij}/2$：其餘埠的反射會再次耦合回 victim。可用網路方程求出接上全部實際阻抗後的 $H_{ij}$。例如網路電流 $\mathbf i$ 取流入網路，外接 Norton 源為 $\mathbf i_N$、外接導納為 $\mathbf Y_{\mathrm{ext}}$，可聯立

$$
\begin{bmatrix}
\mathbf Y_{\mathrm{ext}}&\mathbf I\\
\mathbf I-\mathbf S&-Z_0(\mathbf I+\mathbf S)
\end{bmatrix}
\begin{bmatrix}\mathbf v\\\mathbf i\end{bmatrix}
=\begin{bmatrix}\mathbf i_N\\\mathbf0\end{bmatrix}.
$$

對每個有源 Thévenin 分支，$i_{N,j}=V_{Tj}/R_j$；電容負載須包含在頻率相依 $\mathbf Y_{\mathrm{ext}}$ 內。最後使用有限上升時間、足夠長且有暖機區間的資料，讓 aggressors 有符合預期的碼型與相對延遲，再比較 victim 眼圖。把所有 Tx 設為相同 PRBS、相同相位，可能嚴重低估串擾；不同 seeds 可避免完全重疊，但同一 maximal LFSR 的非零 seeds 本質是同一週期序列的相位平移，並非嚴格統計獨立，仍需檢查相關性及多種相對時序。[L11, p.10–14、19–24]

<a id="q58"></a>

### Q58 Why is FEXT of microstrips usually negative?

先固定符號為 Q52 的端口配置、兩線沿相同方向、victim 電壓同樣對地取正，並討論 aggressor 的**正向上升邊緣**。此時 $v_{\mathrm{FEXT}}\approx(\tau_o-\tau_e)x'/2$；常見表面 microstrip 的 odd 模態較多電場分布在空氣中，有效介電常數通常較低，因而 $v_o>v_e$、$\tau_o<\tau_e$，所以得到負向 FEXT。下降邊緣的極性反過來。這個解釋是模態速度差，不是「有損耗，所以電壓變負」。[L11, p.9、17]

用每單位長度參數也可看出極性。明確定義

$$
\mathbf L=\begin{bmatrix}L_s&L_m\\L_m&L_s\end{bmatrix},\qquad
\mathbf C=\begin{bmatrix}C_s&-C_m\\-C_m&C_s\end{bmatrix},
\quad L_m,C_m>0,
$$

其中 $C_s$ 是 Maxwell 電容矩陣對角元，$-C_m$ 是非對角元，不是把對地電容與互電容的定義混用。令 $k_L=L_m/L_s$、$k_C=C_m/C_s$、$\tau_0=\ell\sqrt{L_sC_s}$，弱耦合時 $\tau_o-\tau_e\approx\tau_0(k_C-k_L)$；因此 $k_L>k_C$ 對應這裡的負向 FEXT。均勻介質下兩模態速度相等，一階容性、感性貢獻可相消；embedded microstrip、非均勻介質、不同幾何及參考極性則可能改變結果。故不能脫離座標、端接與結構，就宣稱「所有 microstrip 的 FEXT 一定為負」。[L11, p.17]

<a id="q59"></a>

### Q59 What is LFSR, its length, seed, and tap polynomials?

LFSR 是 linear feedback shift register，用 $m$ 個儲存位元及 GF(2) 上的 XOR 回授，讓下一狀態成為目前狀態的線性函數。「length $m$」通常指暫存器級數，不是輸出序列的週期；seed 是起始的 $m$ 位元狀態。對採 XOR 的 maximal-length LFSR，非零狀態依序遍歷一次，週期為 $2^m-1$；全零 seed 會永遠卡在零。有限狀態使它必定週期性，輸出是確定的偽隨機序列，不是真正隨機。[L11, p.22–23]

回授可用 GF(2) 多項式描述。若採明確遞迴慣例 $P(x)=x^m+\sum_{i=0}^{m-1}p_ix^i$，則 $u_{k+m}=\bigoplus_{i:p_i=1}u_{k+i}$；例如 primitive polynomial $x^4+x^3+1$ 對應 $u_{k+4}=u_{k+3}\oplus u_k$，使用任一非零四位元 seed 可得週期 15。**primitive** 比 irreducible 更強，才能保證最大週期；任意隨手選的 taps 不保證做到 $2^m-1$。[L11, p.24]

讀 tap 表時必須核對 shift 方向、bit 編號、輸出取哪一級、Fibonacci 或 Galois 結構，以及表列多項式是否採 reciprocal convention；相同多項式名稱不代表可原封不動複製同一組陣列索引。最大長度序列常用於 PRBS 測試及 scrambler，但相同 primitive LFSR 的不同非零 seeds 只是序列相位不同，不能把它們稱作數學上獨立的亂數來源。模擬串擾時還要檢查所用紀錄長度、跨通道相關性與相對時間偏移。[L11, p.11、14、22–24]

<a id="topic-10"></a>

## Jitter basics

<a id="q60"></a>

### Q60 How many different kinds of jitter are there?

沒有一個不依分類方式而變的固定種類數。常見第一層分成 **random jitter（RJ）** 與 **deterministic jitter（DJ）**；再依來源區分 periodic jitter（PJ）、duty-cycle distortion（DCD）、由資料與 ISI 造成的 data-dependent jitter（DDJ），以及與受害資料不相關的 bounded uncorrelated jitter（BUJ，例如部分串音）。這些標籤描述不同面向，並非全部互斥：DCD 可屬於資料相關的 DJ；串音也可能呈現週期性或資料相關性。[L12, pp. 2–9]

所有 timing jitter 都須相對某個參考時脈定義：$J_n=t_n-t_{n,\mathrm{ideal}}$，單位是秒或 UI。固定理想時脈、轉送時脈和經 CDR 恢復的時脈，會得到不同的殘餘 jitter。若數個 jitter 分量可相加且彼此獨立，其 PDF 作摺積、變異數相加；有相關性時須保留 covariance，不能直接使用平方和開根號。[L12, pp. 2, 11–14]

<a id="q61"></a>

### Q61 What is periodic jitter (PJ) and its PDF?

PJ 是隨時間週期變動的邊緣時間誤差，常來自電源切換、鄰近週期訊號或 SSC。以正弦 PJ 為例，$J(t)=A\sin(2\pi f_jt+\theta)$。若觀察涵蓋足夠週期且取樣相位均勻分布，其 PDF 為 **arcsine distribution**：

$$
f_J(x)=\begin{cases}\dfrac{1}{\pi\sqrt{A^2-x^2}},&|x|<A,\\0,&|x|>A.\end{cases}
$$

這裡 $A$ 是時間偏移的峰值，peak-to-peak 為 $2A$、RMS 為 $A/\sqrt2$。兩端密度較高，因為正弦波在極值附近變化較慢；端點發散可積分，不是兩個 Dirac impulses。相同峰值的對稱三角形 PJ 有均勻 PDF，50% 占空比的矩形 PJ 則是兩個等權 impulses。若 jitter 與資料時脈鎖定成有限個取樣相位，實測分布未必等於連續 arcsine PDF。[L12, pp. 5–7]

<a id="q62"></a>

### Q62 What is DCD and its PDF?

DCD（duty-cycle distortion）指高、低邏輯維持時間不對稱，或上升、下降緣相對理想時刻具有不同偏移；可能由時脈占空比、上升／下降時間或判決門檻不對稱造成。將平均偏移扣除後，最簡模型把兩類邊緣放在 $-D/2$ 與 $+D/2$，其中 $D$ 是兩群時間位置的間距。[L12, p. 8]

若兩類邊緣在取樣集合中等權且各自沒有其他 jitter，

$$
f_{\mathrm{DCD}}(x)=\tfrac12\delta(x+D/2)+\tfrac12\delta(x-D/2).
$$

不等權時要改成 $w\delta(x-\mu_1)+(1-w)\delta(x-\mu_2)$；再加入獨立 RJ，便是兩個移位 Gaussian 的混合。注意 DCD 的百分比占空比指標與以 ps 表示的邊緣間距不是同一個定義，換算前要確認使用的是哪一個時脈週期與量測慣例。[L12, pp. 8, 11]

<a id="q63"></a>

### Q63 What is random jitter (RJ) and its PDF?

RJ 是由隨機擾動造成的 timing jitter，常以熱雜訊等許多微小擾動的合成來解釋。在常用的零均值 Gaussian 模型中，$J\sim\mathcal N(0,\sigma^2)$，

$$
f_J(x)=\frac{1}{\sigma\sqrt{2\pi}}\exp\!\left(-\frac{x^2}{2\sigma^2}\right).
$$

$\sigma$ 是 RMS jitter。Gaussian 的支撐範圍無界，因此理論上沒有有限的最壞 peak-to-peak RJ；有限紀錄量到的最大值會隨樣本數改變，應報告 RMS、觀察長度或指定 BER 下的等效 jitter。Gaussian 是模型假設，不能把每一種隨機擾動都當成 Gaussian；而且「每個樣本的邊際分布為 Gaussian」也不代表樣本彼此獨立或頻譜為白色。[L12, pp. 9–10]

<a id="q64"></a>

### Q64 Jitter time function vs jitter PDF

Jitter time function／sequence 保留誤差的**發生順序**：例如 $t_n=nT_b+J(nT_b)$，可直接用來移動模擬中的轉換時刻，也能計算週期、頻譜及自相關。Jitter PDF 則描述「落在某時間偏移附近的機率密度」，適合統計眼圖、組合獨立 jitter 與估算尾端機率。[L12, pp. 4, 10–12]

PDF 不保留時間順序。不同頻率的正弦 PJ，只要振幅相同且相位充分取樣，具有相同 arcsine PDF；把一組 jitter 樣本打亂，直方圖不變，相關性與對 CDR 的影響卻可能改變。因此不能僅由 PDF 還原唯一的 jitter time function，也不能由相同 PDF 推論相同的動態鏈路行為。[L12, p. 10]

<a id="q65"></a>

### Q65 Jitter PDF vs jitter CDF

PDF $f_J(x)$ 是機率密度，單位為 $1/\mathrm{s}$；小區間機率約為 $f_J(x)\,dx$，單一點的密度不是該點的機率。CDF $F_J(x)$ 是無單位的累積機率：

$$
F_J(x)=P(J\le x)=\int_{-\infty}^{x}f_J(u)\,du,
\qquad P(a<J\le b)=F_J(b)-F_J(a).
$$

因此 CDF 非遞減，從 0 趨向 1；右尾機率為 $1-F_J(x)$。PDF 方便看分布形狀，CDF／尾積分才直接對應邊緣提早或延遲超過容許量的機率。從眼圖取得分布時須正規化並區分「有轉換時的條件 jitter 分布」與每個 bit 的轉換機率，不能把未正規化的門檻切片直接當成 CDF。[L13, pp. 2–4, 8–11]

<a id="q66"></a>

### Q66 CDF of Gaussian distribution

若 $J\sim\mathcal N(\mu,\sigma^2)$，令 $z=(x-\mu)/\sigma$，則

$$
F_J(x)=\Phi(z)
=\frac12\left[1+\operatorname{erf}\!\left(\frac{x-\mu}{\sqrt2\sigma}\right)\right]
=\frac12\operatorname{erfc}\!\left(\frac{\mu-x}{\sqrt2\sigma}\right).
$$

右尾用 $Q(z)=1-\Phi(z)=\tfrac12\operatorname{erfc}(z/\sqrt2)$，所以 $P(J>x)=Q((x-\mu)/\sigma)$。檢查符號最簡單的方法是確認 $F_J(\mu)=1/2$，以及 $x\to+\infty$ 時 CDF 趨近 1、右尾趨近 0。計算很小的右尾時直接使用 `erfc` 或 survival function，可避免用 $1-F$ 相減造成的浮點消去誤差。[L13, pp. 3–4]

<a id="q67"></a>

### Q67 What is bit error rate (BER)?

BER 是接收器判決錯一個 bit 的機率；有限紀錄的估計值為

$$
\widehat{\mathrm{BER}}=\frac{N_{\mathrm{wrong}}}{N_{\mathrm{tested}}}.
$$

它取決於整條鏈路及判決方式，包括取樣時刻、門檻、時脈恢復、雜訊與等化器，並非單憑通道的 $S_{21}$ 就能決定。沒有觀察到錯誤只代表此紀錄的計數為零；在獨立同分布 Bernoulli 假設下，測試 $N$ bits 而零錯誤的單側 95% 上界約為 $3/N$，不能宣稱真實 BER 必為零。[L13, pp. 6–7, 14–15]

<a id="q68"></a>

### Q68 How to compute BER from eye?

最直接的方法是在取樣時刻 $t$ 保留**依傳送 bit 分組**的電壓分布。若邏輯 1 的先驗機率為 $p_1$、邏輯 0 為 $p_0$，門檻為 $V_T$，則

$$
\mathrm{BER}(t)=p_1\!\int_{-\infty}^{V_T}f_{V|1,t}(v)\,dv
+p_0\!\int_{V_T}^{\infty}f_{V|0,t}(v)\,dv.
$$

在只考慮邊緣時序誤差、忽略大幅跨越多個 UI 與兩側事件重疊、相鄰邊界錯誤可相加的模型中，對 $t_i<t<t_{i+1}$，講義給出

$$
\mathrm{BER}(t)\simeq\rho_L[1-F_L(t-t_i)]+\rho_R F_R(t-t_{i+1}).
$$

第一項是左邊界太晚，第二項是右邊界太早；分布相同且轉換密度相同時，令 $F_L=F_R=F_J$、$\rho_L=\rho_R=\rho_T$ 即得到講義公式。$\rho_T$ 不一定是 $1/2$。未標記預期 bit 的眼圖直方圖，通常不足以唯一決定 BER；有限眼圖中看不到尾端或眼中心有空白，也不表示極低機率錯誤不存在。抽取 crossing-time histogram 比直接把某一電壓 bin 的密度當成 timing PDF 更明確，因為後者還受 bin 寬與交越斜率影響。[L13, pp. 2, 7–12]

<a id="q69"></a>

### Q69 What information can we get from the bathtub (BER) curve?

Bathtub curve 是在固定門檻及接收器／時脈設定下，掃描取樣相位所得的 $\mathrm{BER}(t)$。兩側靠近轉換邊緣時錯誤率升高，中間通常較低。給定目標 $\mathrm{BER}_0$，左右交點之間的寬度就是該 BER 定義下可用的取樣時間窗，亦可讀出最佳取樣位置、左右 timing margin 及是否有明顯不對稱。[L13, pp. 11–13]

尾端形狀可以協助分析 RJ／DJ，但其解讀依賴模型，不能只靠一條曲線唯一分離所有 jitter 來源。對含理想 Gaussian RJ 的模型，嚴格的「零錯誤最壞 EW」沒有實用的有限尾端界線，因此指定 BER 的 EW 比無機率條件的 peak-to-peak EW 更有意義；每次比較也應使用相同轉換密度、門檻與 CDR 條件。[L13, pp. 5, 13, 16–17]

<a id="q70"></a>

### Q70 How to extrapolate BER curves down to very low probability?

先量測或模擬有足夠事件數的分布／bathtub 尾端，再選擇具有物理依據的模型，例如 bounded DJ 與 Gaussian RJ 的摺積。擬合兩側尾端的平均位置、標準差及權重後，使用模型的 CDF 或 $Q$ 函數計算尚未直接觀察到的低機率區域；Gaussian-dominated 尾端可在合適的 inverse-$Q$／Gaussian-quantile 座標下檢查是否接近直線。[L13, pp. 14–17]

外推是**依模型的估計**，不是新增的實測證據。應交代擬合範圍、事件數、轉換密度、模型殘差與參數不確定性，並檢查其他可觀察尾端區段。稀有 burst、非 Gaussian 雜訊、狀態切換或相關錯誤，都可能使高機率區域的良好擬合失去外推能力；不能把普通的 $\log(\mathrm{BER})$ 圖任意拉直線，也不能用零事件 bin 當成零機率。[L13, pp. 15–18]

<a id="q71"></a>

### Q71 What is the dual-Dirac model and Gaussian mixture model?

Dual-Dirac 把有效 DJ 近似為兩個 impulses，再與獨立 Gaussian RJ 摺積，得到

$$
f_J(x)=\tfrac12\mathcal N(x;\mu_1,\sigma_{RJ}^2)
+\tfrac12\mathcal N(x;\mu_2,\sigma_{RJ}^2).
$$

兩峰間距 $D=\mu_2-\mu_1$ 是尾端擬合的有效 DJ 參數，不必等於所有實際 DJ 的真實 peak-to-peak 總和。常見寫法 $TJ\simeq D+2Q_{\mathrm{target}}\sigma_{RJ}$ 中，$Q_{\mathrm{target}}$ 是對應目標尾機率的 Gaussian 分位倍率，不是可忽略定義的固定常數。例如外側尾由權重 $w$ 的一個 Gaussian 主導時，該尾機率約為 $\epsilon=wQ(q)$，故 $q=Q^{-1}(\epsilon/w)$；若從 BER 換算，還須納入轉換密度與單／雙側誤差的慣例。[L13, p. 16]

更一般的 Gaussian mixture model（GMM）為 $f_J(x)=\sum_{k=1}^{K}w_k\mathcal N(x;\mu_k,\sigma_k^2)$，其中 $w_k\ge0$、$\sum_kw_k=1$，可使用不同平均、變異數與權重擬合非對稱分布。參數更多不代表尾端更可信；有限資料可能無法唯一辨識小權重的大變異數分量，而這種分量可能主導極低 BER 外推，必須檢查可辨識性與模型穩定性。[L13, pp. 17–18]

<a id="topic-11"></a>

## Differential signaling (2)

<a id="q72"></a>

### Q72 Inter-pair skew and intra-pair skew

**Intra-pair skew** 是同一差動對內 P、N 兩路的傳播時間差；它使兩個相反訊號不能同時抵達，導致差動波形失真、DM-to-CM conversion 與眼圖劣化。**Inter-pair skew／lane-to-lane skew** 則是不同差動對或 lanes 之間的到達時間差，影響多 lane 資料對齊，有些接收器可用 deskew buffer 吸收一定範圍。[L14, pp. 2–4]

兩者應以**電氣延遲**衡量，而不只是幾何線長。等長線路仍可能因介電環境、via、連接器及不同層別而有 skew；實際容許值取決於協定、世代、通道類型與接收器能力，不能把某一條 cable 規格的數字直接當成所有 PCB 的規則。[L14, pp. 3–4, 8–10]

<a id="q73"></a>

### Q73 Measurement of skew in practice

時間域方法對兩路施加一致激勵、量取指定電壓門檻的 crossing-time difference。Inter-pair skew 比較兩個差動輸出；intra-pair skew 比較同一差動對的 P、N 單端輸出，須以對應的相反極性門檻或先反相 N 路後比較。講義以 DP 的 15% 與 USB 的 50% crossing 作例子，實作須遵循指定版本的量測定義，並校正探棒／儀器通道自身的延遲。[L14, p. 5]

頻率域可用 VNA 取得差動激勵到兩個單端輸出的傳輸 $H_P=(S_{21}-S_{23})/\sqrt2$、$H_N=(S_{41}-S_{43})/\sqrt2$，在 N 路補償差動訊號既有的 $180^\circ$ 反相後，比較展開相位所對應的延遲。例如

$$
\tau_P(f)-\tau_N(f)=-\frac{\operatorname{unwrap}\arg[H_P(f)/(-H_N(f))]}{2\pi f}.
$$

相位分支與參考面要一致；否則理想對稱線對也會算出假的 skew。講義的 EIPS 再依指定頻帶與權重整合 $\mathrm{skew}(f)$；它是特定的加權相位延遲指標，不能任意換成群延遲導數或自行挑權重。[L14, pp. 6–7]

<a id="q74"></a>

### Q74 What is the fiber-weave effect?

FR-4 類基板由玻璃纖維束與樹脂構成，兩者介電常數不同，形成空間非均勻、常與織紋方向相關的介電環境。若 P 路長時間位於玻璃束上方，N 路較多位於樹脂區，兩條線即使幾何等長，也可能有不同的有效介電常數與傳播速度，產生 intra-pair skew。[L14, pp. 8–10]

以準 TEM 近似表示，$\tau\simeq l\sqrt{\epsilon_{\mathrm{eff}}}/c$，因此 $\Delta\tau\simeq(l/c)(\sqrt{\epsilon_{\mathrm{eff},P}}-\sqrt{\epsilon_{\mathrm{eff},N}})$。當走線寬度小、路徑長且沿著某一纖維束方向持續前進時，局部差異較難平均，可能累積成顯著 skew、模式轉換與差動插入損耗劣化。[L14, pp. 8–10]

<a id="q75"></a>

### Q75 How to mitigate the fiber-weave effect?

主要策略是讓兩條線沿途平均取樣相似的玻璃／樹脂分布：採用較均勻或較密的織法、展纖材料，選擇玻璃與樹脂介電差異較小的材料，或讓走線相對織紋旋轉、採斜向／適當 zig-zag routing。講義列出 10°–20° 旋板與 45° 走線等例子；角度是製程與布局選項，不是普遍保證值。[L14, pp. 11–12]

較寬走線也能平均更多局部介電差異，但仍須維持目標阻抗與空間限制。選擇密織材料、旋轉 panel 或路由方向應與板廠協調；任何補償都需檢查新增線長、彎折、串音及阻抗不連續。單靠尺量等長，無法消除未知的沿線有效介電差異。[L14, pp. 8–12]

<a id="q76"></a>

### Q76 How to compute the mixed-mode S parameters for multiple diff-pairs?

先明確列出每一個 modal port 的單端正、負埠 $(p_i,n_i)$，再把每一對的 power waves 轉為 $a_{c_i}=(a_{p_i}+a_{n_i})/\sqrt2$、$a_{d_i}=(a_{p_i}-a_{n_i})/\sqrt2$，$b$ 同理。對相同實數單端參考阻抗 $Z_0$，將所有 pair 的和／差列放入正交矩陣 $T$，即有

$$
\mathbf a_m=T\mathbf a,\qquad \mathbf b_m=T\mathbf b,
\qquad S_m=TST^{\mathsf T}.
$$

例如有 $N$ 條兩端差動連線時，共有 $4N$ 個單端埠、$2N$ 個成對端點，仍轉成一個 $4N\times4N$ 矩陣。快速計算任一差動傳輸項可用 $S_{dd,ij}=(S_{p_ip_j}-S_{p_in_j}-S_{n_ip_j}+S_{n_in_j})/2$，再依兩端所在位置判斷它是 insertion、NEXT 或 FEXT。輸出模式在第一個下標、輸入模式在第二個下標。[L14, pp. 14–16]

此轉換對應 $Z_d=2Z_0$、$Z_c=Z_0/2$。若測試定義要求其他差動參考阻抗，必須先作一致的 renormalization；例如講義的 85-ohm DM 作法是先將 SE 改為 42.5 ohm，再轉成 DM 85 ohm／CM 21.25 ohm。遇到不同單端參考阻抗或另有指定的 CM 阻抗時，不能不加調整就沿用這個簡單正交式。[L14, pp. 16–17]

<a id="q77"></a>

### Q77 How is crosstalk gauged in USB/PCIe standards?

共同概念是先使用正確埠配對與參考阻抗的 **differential crosstalk S parameters**，再依指定頻帶、激勵／接收端權重及 aggressor 集合評估。講義引用的 USB Type-C Cable and Connector Specification release 2.2（2022）以 integrated NEXT／FEXT 為例，其形式是

$$
I_{XT}=10\log_{10}\!\left[
\frac{\int_0^{f_{\max}}W(f)|S_{dd}^{XT}(f)|^2\,df}
{\int_0^{f_{\max}}W(f)\,df}\right],\quad XT\in\{NEXT,FEXT\},
$$

其中權重與規定的輸入脈衝頻譜有關。IFEXT 必須使用 FEXT；講義該列仍寫 NEXT 應視為排版錯字。講義引用的 PCIe 5.0／6.0 internal-cable 草案則示範 power-sum NEXT／FEXT，例如 $PSNEXT_i=10\log_{10}\sum_j|S_{dd,ij}^{NEXT}|^2$，配合頻率 mask 與 ccICN 加權積分雜訊判準。ccICN 的結果是由加權積分得到的 RMS 串音電壓，與前述正規化 integrated-crosstalk 的 dB 比值不同。功率相加不同於把所有 aggressor 的複數電壓同相相加；後者還涉及彼此相關性。[L14, pp. 18–19; W03]

這些是**特定版本與通道類型**的例子，不是所有 USB／PCIe 世代共用的單一限制。實際合規須逐項核對指定版本、測試件、參考阻抗、權重、頻帶、aggressor 組合及 passing criteria；不能只檢查某個頻點的串音，也不能把講義示例的 dB 數字當成目前所有產品的門檻。[L14, pp. 17–19]

<a id="q78"></a>

### Q78 Techniques to reduce differential crosstalk?

先降低不必要的電磁耦合：增加相鄰線對間距、縮短平行共走距離、讓參考平面更接近、保持連續回流路徑，並在相鄰訊號層適當採用正交走線。維持每對 P/N 的對稱性與成對路由，有助於減少外洩場及共模轉換；放慢邊緣則須以仍滿足資料率與眼圖為前提。均勻介質 stripline 可降低典型的 FEXT 失配，但不能保證完整實體通道的差動串音恆為零。[L14, pp. 20–22]

Guard trace 只有在良好接地、適當 via stitching 且沒有形成長浮接諧振段時才有幫助。更進階的方法利用不同耦合路徑的極性抵消：switched／jogged routing、edge 與 broadside coupling 的組合、最佳線對／via 偏移等。這類抵消受幅度、相位及頻率影響，須用完整 EM S 參數及多 aggressor 時域／眼圖驗證，不能只依單一頻點或幾何直覺判定成功。[L14, pp. 21, 23–31]

<a id="topic-12"></a>

## PAM4 signaling

<a id="q79"></a>

### Q79 What is PAM4 signaling and why people use it?

PAM4 使用四個振幅階，每個 symbol 理想上攜帶 $\log_2 4=2$ bits；典型對稱電壓為 $\{-1,-1/3,1/3,1\}V_{pp}/2$。NRZ／PAM2 每個 symbol 只有 1 bit，所以在**相同線路 bit rate** 下，PAM4 的 baud rate 為 NRZ 的一半、symbol interval 為兩倍，所需的主要頻寬也較小。[L15, pp. 2–4]

這讓高頻損耗大的通道可在較低 baud rate 傳送更多資料。代價是固定外側 $V_{pp}$ 時，相鄰階距只有 $V_{pp}/3$，需要三個判決門檻，對雜訊、非線性及電壓誤差更敏感。若比較的是**相同有效資料率**，還須計入編碼與 FEC 效率 $\eta$：$R_{\mathrm{net}}=R_s\log_2(M)\eta$；只有效率相同時，才可直接說 baud rate 恰好減半。[L15, pp. 3–5, 11–12]

<a id="q80"></a>

### Q80 PSD of PAM4 signals?

對獨立、等機率、零均值的 PAM4 symbols，若幅度為 $a_n\in\{-1,-1/3,1/3,1\}V_{pp}/2$，則 $\sigma_a^2=E[a_n^2]=5V_{pp}^2/36$。令單位高度的 pulse 為 $p(t)$，$x(t)=\sum_na_np(t-nT_s)$，time-averaged、雙邊電壓 PSD 為

$$
S_x(f)=\frac{\sigma_a^2}{T_s}|P(f)|^2
=\frac{5V_{pp}^2T_s}{36}\operatorname{sinc}^2(fT_s)\operatorname{sinc}^2(fT_r),
$$

最後一式使用與講義相同的有限線性 rise-time pulse，$\operatorname{sinc}(x)=\sin(\pi x)/(\pi x)$。此式單位為 $\mathrm{V^2/Hz}$；對電阻負載的功率 PSD 要再除以相應電阻。非零均值、非獨立或非等機率符號需加入均值／相關性的影響，不能直接使用這個簡式。[L15, pp. 7–8]

在**相同 baud rate、相同 $V_{pp}$ 與 pulse shape** 下，PAM4 PSD 是 NRZ 的 $5/9$；若改成相同平均功率正規化，這個幅度比例也會改變。比較相同 bit rate 時，$T_{s,PAM4}=2T_{b,NRZ}$，主要頻譜壓縮至較低頻率，不能在每個頻點都套用固定 $5/9$。講義 p. 8 的 NRZ 簡寫應維持兩個 sinc 的平方，與 p. 2 及一般 PAM PSD 公式一致。[L15, pp. 2, 7–8]

<a id="q81"></a>

### Q81 Eye diagram of PAM signals?

理想 $M$-PAM 有 $M$ 個電壓階、$M-1$ 個相鄰判決區，因此 PAM4 有三個眼。每個眼都需要分別檢查 EH、EW、jitter、門檻及目標錯誤率，不能只展示最好的中間眼。在等機率、相鄰階具有相同 Gaussian 雜訊的理想情況，最佳門檻位於相鄰電平中點；電平間距、雜訊或機率不對稱時，門檻未必是等分的。[L15, pp. 4–6, 11]

固定外側 $V_{pp}$ 的理想 PAM4 相鄰階距是 NRZ 的 $1/3$，對應電壓間距比 $20\log_{10}3\approx9.54$ dB。PAM4 有更多不同跨度的轉換；若各轉換都使用相同的有限 rise/fall duration，不同跨度穿越同一門檻的時刻可能不同，造成資料相關的 crossing-time spread。這項展寬取決於轉換波形，不是所有理想化 PAM4 波形都必然具有相同的「固有 jitter」。[L15, pp. 5–6]

<a id="q82"></a>

### Q82 When will PAM4 be better than NRZ?

當通道的高頻損耗非常嚴重，降低 baud rate 所恢復的訊號幅度與時間裕度，可能超過 PAM4 的較小階距及額外接收器代價。比較相同線路 bit rate $R_b$ 時，NRZ 與 PAM4 的 Nyquist frequencies 分別為 $R_b/2$ 與 $R_b/4$。若以正值衰減量 $IL(f)=-20\log_{10}|S_{dd21}(f)|$ 表示，講義提出經驗觀察：$IL(R_b/2)-IL(R_b/4)$ 超過約 9.5 dB 時，PAM4 的眼高可能較有利，與固定 $V_{pp}$ 下 $1/3$ 的階距代價相呼應。[L15, pp. 9–10]

這是所研究通道與激勵條件下的**經驗規則**，不是只看兩個 S-parameter 數值就能判定 BER 的定理。實際選擇須在相同有效資料率、功率／擺幅限制、編碼與 FEC 條件下，比較三個眼、雜訊增強、ISI、線性度、jitter、等化能力及接收器複雜度；低損耗、雜訊主導的通道可能仍較適合 NRZ。[L15, pp. 5, 9–12]

<a id="topic-13"></a>

## Spectrum control by scrambler and SSC

<a id="q83"></a>

### Q83 Characterization of signal spectra

週期訊號可用 Fourier series 描述，功率集中在離散 harmonics；平穩隨機訊號則以 PSD 描述平均功率如何分布於頻率。實際訊號可同時含有離散 tones 與連續背景，例如重複資料模式加上雜訊；因此「random」不自動代表沒有 DC 或任何譜線，非零均值與週期／相關結構都可能留下離散成分。[L16, pp. 2–3]

量測 EMI 時應區分**一條譜線的總功率**與 **每 Hz 的 PSD**。若某一未展頻 tone 完整落在接收器帶寬內，其積分功率不會隨 RBW 等比例改變；平坦寬頻 PSD 在量測帶寬內的功率則近似 PSD 乘上等效雜訊帶寬。集中譜線常造成較高的窄頻峰值，但輻射／干擾仍取決於耦合與接收器頻帶；把能量攤寬不代表所有 RFI 情況都改善。[L16, pp. 2–3, 20]

<a id="q84"></a>

### Q84 What is scrambler and why people use it?

Scrambler 是可逆的資料擾亂器，用來降低常見資料中的重複結構、偏壓或不利轉換特性，讓傳送資料在統計上較接近隨機序列。Additive scrambler 將資料與 LFSR 產生的 PRBS 作 XOR：$s[n]=x[n]\oplus p[n]$；接收端用相同且同步的 PRBS 再 XOR 還原，因此需要一致的多項式、狀態與同步方式。[L16, pp. 4–7]

Self-synchronizing／multiplicative scrambler 則使用先前的輸出位元，例如 $s[n]=x[n]\oplus s[n-18]\oplus s[n-23]$；接收端以收到的序列作同樣的延遲 XOR，經有限歷史填滿後可恢復資料，不需另外對齊一個獨立 PRBS 狀態，但通道錯誤可能擴散成數個輸出錯誤。這是 $GF(2)$ 上的線性運算，不是一般實數電壓的加法，也不是加密。[L16, pp. 8–9]

<a id="q85"></a>

### Q85 Effect of scrambler on signal spectrum?

Scrambling 通常降低重複資料所造成的突出譜線，把能量分散在更多頻率成分中；頻譜包絡仍受 symbol interval、pulse shape 與 rise time 控制。對長週期 PRBS，線與線的間距可能非常小，量測上接近寬頻背景，但有限狀態的 scrambler 與重複輸入仍可能產生週期輸出，不能把「看起來較隨機」說成數學上真正非週期。[L16, pp. 3–4, 10]

Scrambler **不保證所有輸入都被白化、DC 完全消失或具有固定的最大連續相同位元數**。例如 additive scrambler 的輸入若恰好等於其 PRBS，XOR 輸出就是全零。要保證特定 run-length 或 DC 限制，必須依賴另有約束的編碼／協定設計；scrambling 本身主要提供統計上的改善。[L16, pp. 4, 8–10]

<a id="q86"></a>

### Q86 What is spread-spectrum clocking (SSC)?

SSC 讓時脈的瞬時頻率在一個小範圍內緩慢變動，以 frequency modulation 將原本集中的 clock harmonics 分散到較寬的頻帶，降低指定量測帶寬看到的峰值。它保留可用的時脈功能；與資料 scrambling 不同，SSC 調整的是邊緣時間與時脈頻率，而不是將時脈邏輯狀態隨機 XOR。[L16, pp. 11–14]

SSC 通常沒有等比例降低訊號總功率，而是重新分配頻譜。代價是增加可預期的頻率偏移與 timing modulation，因此須符合時脈容差、PLL／CDR 追蹤能力及協定限制；窄頻 EMI 峰值改善，也不保證鄰近無線接收器在所有條件下都受益。[L16, pp. 11–14, 20]

<a id="q87"></a>

### Q87 Spreading rate, modulation rate, spreading waveform, spreading style

令 $f_c$ 為 nominal clock frequency，$\Delta f=f_{\max}-f_{\min}$ 為**完整展頻寬度**，則 fractional spreading range 可寫成 $\delta_{pp}=\Delta f/f_c$。Modulation rate 是頻率變化一循環的速率 $f_m=1/T_m$；modulation waveform 描述 $f(t)$ 如何隨時間改變，例如三角形、正弦或經修形的 profile。[L16, pp. 12–14]

若使用上述同一個完整寬度定義，down-spread 範圍為 $[f_c-\Delta f,f_c]$，center-spread 為 $[f_c-\Delta f/2,f_c+\Delta f/2]$，up-spread 為 $[f_c,f_c+\Delta f]$。資料表可能將 center-spread 寫成 $\pm\delta$，此時 $\delta$ 表示半寬；比較百分比或 ppm 前須確認慣例。不同 profile 在各頻率停留的時間不同，會影響展頻均勻度；正弦調變通常較偏重兩端。[L16, pp. 12–13; W02]

<a id="q88"></a>

### Q88 Estimation of EMI reduction by SSC

假設一條原本集中的譜線功率 $P$ 近似均勻地展開在 $\Delta f$，接收器等效帶寬為 $B$，且 $B\ll\Delta f$，則帶內平均功率約為 $P B/\Delta f$，所以峰值降低的初步估算為

$$
R_{\mathrm{EMI}}\simeq10\log_{10}\!\left(\frac{\Delta f}{B}\right)\ \mathrm{dB}.
$$

例如 $\Delta f=3.25$ MHz、$B=100$ kHz，得到約 $15.1$ dB。這是功率／頻寬比的 $10\log_{10}$，不是 $20\log_{10}$。實際應使用量測濾波器相應的等效帶寬；頻譜不平坦、ripple、偵測器與掃描／觀察設定，都會使實際最高峰不同於平均估算。[L16, pp. 15–16; W02]

此估算需充分觀察調變循環，且帶內可把展頻視為近似連續的功率分布。若週期 SSC 的離散 sidebands 已被很窄的 RBW 分開，繼續縮小 RBW 不會使單一譜線功率無限降低；若 $B$ 已涵蓋整個展頻寬度，也沒有上述分攤增益。廠商近似亦要求時脈頻率遠高於調變頻率，量測掃描足夠慢；不能忽略指定的 peak／average／quasi-peak detector 條件而保證同一個抑制數字。[L16, pp. 19–20; W02]

<a id="q89"></a>

### Q89 Simulation of SSC

先指定符合 spread style 的瞬時頻率 $f(t)$，再積分得到相位：

$$
\theta(t)=\theta(0)+2\pi\int_0^t f(u)\,du,
\qquad x(t)=A\cos\theta(t).
$$

數位 clock 可由相位的週期位置產生高低電平，或求出指定相位 crossing 的 edge times；50% duty clock 的相鄰跳變間相位增加 $\pi$。當調變很慢時，可用 $t_{k+1}\approx t_k+1/[2f(t_k)]$ 逐邊緣近似，並檢查時間步長與相位積分誤差。**一般情況不能寫成 $\cos[2\pi f(t)t]$**，因為其瞬時頻率是 $f(t)+t f'(t)$；講義的簡寫必須以真正的相位積分或適用近似解讀。[L16, pp. 15, 17–18]

接著加入所需 rise/fall time，使用足夠長的紀錄算 Fourier／PSD，再以指定 RBW 的濾波器或頻帶積分比較展頻前後。若要用一個 $T_m$ 的 Fourier series，除了 $f(t)$ 週期性，還必須確認相位在一循環內累積整數個 $2\pi$、波形確實接續；否則應用較長紀錄與適當 window。FFT bin spacing $1/T_{obs}$ 也不自動等於真實儀器 RBW。最後檢查平均頻率、頻率上下限、調變週期、功率守恆與峰值降低是否符合設定。[L16, pp. 17–20]

<a id="lecture-sources"></a>

## 講義對照

課程講義、原始回顧題與教學材料來源致謝：**國立中央大學電機工程學系周求致教授（Prof. Chiu-Chih Chou, Taiwan）**。詳見[來源署名與連結](../../ATTRIBUTION.md)。本頁逐題解答為另行整理的複習內容。

下列代號對應課程提供的講義檔案。檔名保留原拼字，頁碼均為 PDF 頁次。

| 代號 | 講義檔名 | 主要對應題目 |
| --- | --- | --- |
| L01 | `1 Time-freq transform.pdf` | Q01–Q07 |
| L02 | `2 Time-freq transform.pdf` | Q08–Q12 |
| L03 | `3 Time-freq transform.pdf` | Q13–Q15 |
| L04 | `4 Spectrum of digital signal.pdf` | Q16–Q22 |
| L05 | `5 Eye diagram.pdf` | Q23–Q26 |
| L06 | `6 Superposition ISI PDA.pdf` | Q27–Q31 |
| L07 | `7 CTLE.pdf` | Q32–Q35 |
| L08 | `8 FFE.pdf` | Q36–Q40 |
| L09 | `9 DFE.pdf` | Q41–Q43 |
| L10 | `10 Differential signaling (1).pdf` | Q44–Q51 |
| L11 | `11 Crosstalk (1).pdf` | Q52–Q59 |
| L12 | `12 Jitter basics.pdf` | Q60–Q64 |
| L13 | `13 Jitter basics (2).pdf` | Q65–Q71 |
| L14 | `14 Differential siging and crosstalk (2).pdf` | Q72–Q78 |
| L15 | `15 PAM4.pdf` | Q79–Q82 |
| L16 | `16 Scrambler & SSC.pdf` | Q83–Q89 |

## 補充參考

- **W01**：[Keysight，About Resolution Bandwidth](https://helpfiles.keysight.com/csg/89600B/Webhelp/Subsystems/gui/content/rbw_about_rbw.htm)。用於區分 RBW 的定義與等效雜訊頻寬；不同分析儀的 RBW 標示不一定相同。
- **W02**：[Analog Devices，Clock Generation with Spread Spectrum](https://www.analog.com/en/resources/technical-articles/clock-generation-with-spread-spectrum.html)。用於 SSC 的展頻範圍、調變率與頻譜峰值降低估算。
- **W03**：[Ansys，Calculating Integrated Crosstalk Noise](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Circuit/Content/Circuit/CalculatingIntegratedCrosstalkNoise.htm)。用於說明 ICN 方法與 NEXT/FEXT aggressor 設定；實際測試條件須依指定介面與版本。

USB 與 PCIe 的數值門檻須查使用中的規格版本與測試文件；可從 [USB Type-C 規格頁](https://www.usb.org/usb-type-cr-cable-and-connector-specification)與 [PCI-SIG 規格目錄](https://pcisig.com/specifications)查找。本回顧說明講義中的評估方法，並未把某個版本的門檻當成所有世代的共同要求。
