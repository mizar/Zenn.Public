---
title: "除算高速化のための前計算"
emoji: "🐙"
type: "tech"
topics: []
published: false
---


## 除算の商を求める3つのアプローチ

この節では、除算の商 $\lfloor x/d\rfloor$ を求めるアプローチを3つ考えてみます。

- a. 丁度の値を求める
    - $\lfloor x/d\rfloor=f_d(x)$
- b. 丁度か $1$ 小さい値を求める
    - $\lfloor x/d\rfloor-1\leq f_d(x)\leq\lfloor x/d\rfloor$
- c. 丁度か $1$ 大きい値を求める
    - $\lfloor x/d\rfloor\leq f_d(x)\leq\lfloor x/d\rfloor+1$

コンパイラによる最適化で使われるパターンは主に 「a. 丁度の値を求める」 ですが、

後者の2つのアプローチは、求めた除算の商の近似値から、仮に除算の剰余を計算してみて、商が大きすぎたり小さすぎたりしたと分かれば商の近似値から1だけ調整して、正しい除算の商とするものです。後ほど、 「b. 丁度か $1$ 小さい値を求める」 パターンで 被除数 $0\leq x\lt 2^{128}$ に対して 除数 $d=10^{32}$ での除算最適化を行う例を紹介します。

また、以下のような略記を所々で使うので定義しておきます。

- `mulh(m, x)` $\displaystyle=\left\lfloor\frac{mx}{2^W}\right\rfloor\quad\left(m,x\in\mathbb{Z},\quad 0\le m\lt 2^W,\quad 0\le x\lt 2^W\right)$
    - $W$ ビット整数どうしの乗算の積の、上位 $W$ ビットを取り出す関数。（主に $W=64$）
- $\displaystyle\left(-2^t\right)\operatorname{mod}d=(d-1)-\left(\left(2^t-1\right)\operatorname{mod}d\right)\quad\left(d,t\in\mathbb{Z},\quad 1\leq d\lt 2^W,\quad 0\leq t\right)$
    - 負の被除数に対する剰余計算。ここでは主に、「$((2^t\text{ 以上で }d\text{ の倍数となる最小の整数}) - 2^t)$ はいくつ?」という文脈で使います。
    - 逆に、「$(2^t-(2^t\text{ 以下で }d\text{ の倍数となる最大の整数}))$ はいくつ?」という文脈では、 $2^t\operatorname{mod}d$ のように表記します。

### a. 丁度の値を求める ($d:$ 除数, $x:$ 被除数, $W:$ 整数型のbit長)

$\displaystyle d,r,t,x\in\mathbb{Z},\quad 1\leq d,\quad 0\leq t,\quad 0\leq r=\left(-2^t\right)\operatorname{mod}d\lt d,\quad 0\le rx\lt 2^t$

$\displaystyle\quad\longrightarrow\quad\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\times\frac{x}{2^t}\right\rfloor=\left\lfloor\frac{2^t+r}{d}\times\frac{x}{2^t}\right\rfloor=\left\lfloor\frac{x}{d}+\frac{rx}{2^t d}\right\rfloor=\left\lfloor\frac{x}{d}\right\rfloor$

$2^t$ 以上で最小の $d$ の倍数を $(2^t+r)$ として、仮に $0\leq rx\lt 2^t$ だったとすると、 $\displaystyle 0\leq\frac{rx}{2^t}\lt 1,\quad 0\leq\frac{rx}{2^td}\lt\frac{1}{d}$ が導かれることから、 $\displaystyle\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\times\frac{x}{2^t}\right\rfloor=\left\lfloor\frac{x}{d}\right\rfloor$ となります。

#### $2^W$ 未満の整数を扱う場合の例

$W,d,x\in\mathbb{Z},\quad 2\leq W,\quad 1\leq d\lt 2^W,\quad 0\leq x\lt 2^W$ において 除算の商 $\left\lfloor x/d\right\rfloor$ を計算:

- $d=1$ の場合:
    - $\displaystyle\left\lfloor\frac{x}{d}\right\rfloor=x$
    - 擬似コード: `x`
- $\log_2d=\lceil\log_2d\rceil$ ($d$ が $2$ の累乗数) の場合:
    - $\displaystyle t=\log_2d,\quad\left\lfloor\frac{x}{d}\right\rfloor=\left\lfloor\frac{x}{2^{\log_2d}}\right\rfloor$
    - 前計算: $s=\log_2d$
    - 擬似コード: `x >> s`
- $\displaystyle 2\leq d\lt 2^W,\quad t=\lfloor\log_2(d-1)\rfloor+W,\quad\left((-2^t)\operatorname{mod}d\right)x\lt 2^t$ の場合:
    - $0\leq x\lt 2^{W-1}\quad\longrightarrow\quad\left((-2^t)\operatorname{mod}d\right)x\lt 2^t$
    - $0\leq x\lt 2^W,\quad (-2^t)\operatorname{mod}d\lt 2^{(t-W)}\quad\longrightarrow\quad \left((-2^t)\operatorname{mod}d\right)x\lt 2^t$
    - $\displaystyle 2^{(t-W)}\leq(d-1)\lt d\leq 2^{(t-W+1)},\quad 2^{(W-1)}\leq\lceil 2^t/d\rceil\lt 2^W,\quad 0\leq t-W=\lfloor\log_2(d-1)\rfloor\lt W,$
    - $\displaystyle \left\lfloor\frac{x}{d}\right\rfloor=\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\frac{x}{2^t}\right\rfloor=\left\lfloor\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\frac{x}{2^W}\right\rfloor\frac{1}{2^{(t-W)}}\right\rfloor$
    - 前計算: $\displaystyle s=\lfloor\log_2(d-1)\rfloor,\quad t=s+W,\quad m=\left\lceil 2^t/d\right\rceil$
    - 擬似コード: `mulh(m, x) >> s`
- $\left(2\leq d\lt 2^W,\quad 0\leq x\lt 2^W\right)$ の範囲にて、主に上記以外の場合: 
    - $\displaystyle t=\lfloor\log_2(d-1)\rfloor+W+1,\quad$
    - $\displaystyle 2^{(t-W-1)}\leq(d-1)\lt d\leq 2^{(t-W)},\quad 2^W\leq\lceil 2^t/d\rceil\lt 2^{(W+1)},\quad$
    - $\displaystyle y=\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\frac{x}{2^W}\right\rfloor-x=\left\lfloor\left(\left\lceil\frac{2^t}{d}\right\rceil-2^W\right)\frac{x}{2^W}\right\rfloor,\quad 0\leq y\leq x\lt 2^W,\quad$
    - $\displaystyle\left\lfloor\frac{x}{d}\right\rfloor=\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\frac{x}{2^t}\right\rfloor=\left\lfloor\frac{x+y}{2^{(t-W)}}\right\rfloor=\left\lfloor\left(\left\lfloor\frac{x-y}{2}\right\rfloor+y\right)\frac{1}{2^{(t-W-1)}}\right\rfloor$
    - 前計算: $\displaystyle s=\lfloor\log_2(d-1)\rfloor,\quad t=s+W+1,\quad m=\left\lceil 2^t/d\right\rceil-2^W$
    - 擬似コード: `y = mulh(m, x); (((x - y) >> 1) + y) >> s`

#### 例

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{67}\lt\frac{2^{93}}{39898175}\quad\longrightarrow\quad \left\lfloor\left\lceil\frac{2^{93}}{998244353}\right\rceil\times\frac{x}{2^{93}}\right\rfloor=\left\lfloor\frac{2^{93}+39898175}{998244353}\times\frac{x}{2^{93}}\right\rfloor=\left\lfloor\frac{x}{998244353}\right\rfloor$

$\displaystyle x\in\mathbb{Z},\quad 0\le x\lt 2^{70}\lt\frac{2^{90}}{875776}\quad\longrightarrow\quad\left\lfloor\left\lceil\frac{2^{90}}{10^8}\right\rceil\times\frac{x}{2^{90}}\right\rfloor=\left\lfloor\frac{2^{90}+875776}{10^8}\times\frac{x}{2^{90}}\right\rfloor=\left\lfloor\frac{x}{10^8}+\frac{875776 x}{2^{90}\times 10^8}\right\rfloor=\left\lfloor\frac{x}{10^8}\right\rfloor$

$\displaystyle x\in\mathbb{Z},\quad 0\le x\lt 10^{10}\lt 2^{34}\lt\frac{2^{45}}{1168}\quad\longrightarrow\quad\left\lfloor\left\lceil\frac{2^{45}}{10^4}\right\rceil\times\frac{x}{2^{45}}\right\rfloor=\left\lfloor\frac{2^{45}+1168}{10^4}\times\frac{x}{2^{45}}\right\rfloor=\left\lfloor\frac{x}{10^4}+\frac{1168 x}{2^{45}\times 10^4}\right\rfloor=\left\lfloor\frac{x}{10^4}\right\rfloor$

$\displaystyle x\in\mathbb{Z},\quad 0\le x\lt 43690\lt\frac{2^{19}}{12}\quad\longrightarrow\quad\left\lfloor\left\lceil\frac{2^{19}}{10^2}\right\rceil\times\frac{x}{2^{19}}\right\rfloor=\left\lfloor\frac{2^{19}+12}{10^2}\times\frac{x}{2^{19}}\right\rfloor=\left\lfloor\frac{x}{10^2}+\frac{12 x}{2^{19}\times 10^2}\right\rfloor=\left\lfloor\frac{x}{10^2}\right\rfloor$

$\displaystyle x\in\mathbb{Z},\quad 0\le x\lt 1024=\frac{2^{11}}{2}\quad\longrightarrow\quad \left\lfloor\left\lceil\frac{2^{11}}{10}\right\rceil\times\frac{x}{2^{11}}\right\rfloor=\left\lfloor\frac{2^{11}+2}{10}\times\frac{x}{2^{11}}\right\rfloor=\left\lfloor\frac{x}{10}+\frac{2 x}{2^{11}\times 10}\right\rfloor=\left\lfloor\frac{x}{10}\right\rfloor$

$\displaystyle x\in\mathbb{Z},\quad 0\le x\lt 2^{64}\lt\frac{2^{67}}{5}\quad\longrightarrow\quad \left\lfloor\left\lceil\frac{2^{67}}{7}\right\rceil\times\frac{x}{2^{67}}\right\rfloor=\left\lfloor\frac{2^{67}+5}{7}\times\frac{x}{2^{67}}\right\rfloor=\left\lfloor\frac{x}{7}+\frac{5 x}{2^{67}\times 7}\right\rfloor=\left\lfloor\frac{x}{7}\right\rfloor$ $\quad\displaystyle\left(2^{64}\lt\frac{2^{67}+5}{7}\lt 2^{65}\text{ である事に実装上注意}\right)$

### b. 丁度か $1$ 小さい値を求める ($d:$ 除数, $x:$ 被除数, $s,t:$ パラメータ [ $x$ の範囲に影響 ])

$\displaystyle d,s,t,x\in\mathbb{Z},\quad 0\le \min(s,t,x),\quad 2^s\le d\le 2^{(s+t)},\quad\left(2^{(s+t)}\operatorname{mod}d\right)\left(x-2^s+1\right)\lt 2^{(s+t)}\left(d-2^s+1\right)$

$\displaystyle\quad\longrightarrow\quad\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)\le\left\lfloor\left\lfloor\frac{x}{2^s}\right\rfloor\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\frac{1}{2^t}\right\rfloor\le\left\lfloor\frac{x}{d}\right\rfloor$


- 前計算: $\displaystyle m=\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor$
- 擬似コード: `((x >> s) * m) >> t`
- $s=0$ の場合の擬似コード: `(m * x) >> t`
- $s=0,t=W$ の場合の擬似コード: `mulh(m, x)` ※ $d=1$ の場合、 $m=\lfloor 2^W/d\rfloor=2^W$ となる事に注意

#### 導出

$\displaystyle d,s,t,x\in\mathbb{Z},\quad 0\le \min(s,t,x),\quad 2^s\le d\le 2^{(s+t)},\quad\left(2^{(s+t)}\operatorname{mod}d\right)(x-2^s+1)\lt 2^{(s+t)}(d-2^s+1)$

$\displaystyle\quad\longrightarrow\quad 0\leq\left(2^{(s+t)}\operatorname{mod}d\right)x\lt 2^{(s+t)}d-\left(2^s-1\right)\left(2^{(s+t)}-\left(2^{(s+t)}\operatorname{mod}d\right)\right)$

$\displaystyle\quad\longrightarrow\quad\left(\frac{\left(2^s-1\right)\left(2^{(s+t)}-\left(2^{(s+t)}\operatorname{mod}d\right)\right)}{2^{(s+t)}d}-1\right)\lt\left(-\frac{2^{(s+t)}\operatorname{mod}d}{2^{(s+t)}d}x\right)\leq 0 \quad\quad$ $\displaystyle\left(-\frac{1}{2^{(s+t)}d}\text{ を掛ける}\right)$

$\displaystyle\quad\longrightarrow\quad\left(\frac{x}{d}-1\right)\lt\left(\frac{x-2^s+1}{2^s}\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\frac{1}{2^t}\right)=\frac{\left(x-2^s+1\right)\left(2^{(s+t)}-\left(2^{(s+t)}\operatorname{mod}d\right)\right)}{2^{(s+t)}d}\leq\left(\frac{x}{d}-\frac{2^s-1}{2^{(s+t)}}\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\right)\leq\frac{x}{d}\quad\quad$ $\displaystyle\left(\left(\frac{x}{d}-\frac{2^s-1}{2^{(s+t)}}\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\right)=\frac{2^{(s+t)}x-(2^s-1)\left(2^{(s+t)}-\left(2^{(s+t)}\operatorname{mod}d\right)\right)}{2^{(s+t)}d}\text{ を足す},\ \frac{2^{(s+t)}-\left(2^{(s+t)}\operatorname{mod}d\right)}{d}=\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\right)$

$\displaystyle\quad\longrightarrow\quad \left(\frac{x}{d}-1\right)\lt\left(\frac{x-2^s+1}{2^s}\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\frac{1}{2^t}\right)\le\left(\left\lfloor\frac{x}{2^s}\right\rfloor\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\frac{1}{2^t}\right)\le\frac{x}{d}\quad\quad$ $\displaystyle\left(\frac{x}{2^s}-1\lt\frac{x-2^s+1}{2^s}\le\left\lfloor\frac{x}{2^s}\right\rfloor\le\frac{x}{2^s}\right)$

$\displaystyle\quad\longrightarrow\quad\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)\le\left\lfloor\left\lfloor\frac{x}{2^s}\right\rfloor\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\frac{1}{2^t}\right\rfloor\le\left\lfloor\frac{x}{d}\right\rfloor\quad\quad$ $\left(\text{床関数(floor)を適用}\right)$


#### 例

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{64}\lt\frac{2^{64}\times 998244353}{2^{64}\operatorname{mod}998244353}\quad\longrightarrow\quad\left\lfloor\frac{x}{998244353}\right\rfloor-1\leq\left\lfloor\left\lfloor\frac{2^{64}}{998244353}\right\rfloor\times\frac{x}{2^{64}}\right\rfloor=\left\lfloor\frac{2^{64}-932051910}{998244353}\times\frac{x}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{998244353}\right\rfloor,\quad \left\lfloor\frac{2^{64}}{998244353}\right\rfloor=18479187002$

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{128}\displaystyle\quad\longrightarrow\quad\left(\left\lfloor\frac{x}{10^{32}}\right\rfloor-1\right)\leq\left\lfloor\left\lfloor\frac{x}{2^{64}}\right\rfloor\times\left\lfloor\frac{2^{128}}{10^{32}}\right\rfloor\times\frac{1}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{10^{32}}\right\rfloor,\quad\left\lfloor\frac{2^{128}}{10^{32}}\right\rfloor=3402823$

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 10^{32}\displaystyle\quad\longrightarrow\quad\begin{cases}\displaystyle\left(\left\lfloor\frac{x}{10^{16}}\right\rfloor-1\right)\leq\left\lfloor\left\lfloor\frac{x}{2^{43}}\right\rfloor\times\left\lfloor\frac{2^{107}}{10^{16}}\right\rfloor\times\frac{1}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{10^{16}}\right\rfloor\\\displaystyle\left(\left\lfloor\frac{x}{10^{16}}\right\rfloor-1\right)\leq\left\lfloor\left\lfloor\frac{x}{2^{48}}\right\rfloor\times\left\lfloor\frac{2^{112}}{10^{16}}\right\rfloor\times\frac{1}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{10^{16}}\right\rfloor\\\displaystyle\left(\left\lfloor\frac{x}{10^{16}}\right\rfloor-1\right)\leq\left\lfloor\left\lfloor\frac{x}{2^{52}}\right\rfloor\times\left\lfloor\frac{2^{116}}{10^{16}}\right\rfloor\times\frac{1}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{10^{16}}\right\rfloor\end{cases}$

「a. 丁度の値を求める」 アプローチで $d=10^{32}, 0\leq x\lt 2^{128}$ について目指す例:

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{128}\displaystyle\quad\longrightarrow\left\lfloor\left\lceil\frac{2^{235}}{10^{32}}\right\rceil\times\frac{x}{2^{235}}\right\rfloor=\left\lfloor\frac{x}{10^{32}}\right\rfloor,\quad 2^{128}\lt\left\lceil\frac{2^{235}}{10^{32}}\right\rceil\lt 2^{129}$

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{128}\displaystyle\quad\longrightarrow\left\lfloor\left\lfloor\frac{x}{2^{32}}\right\rfloor\times\left\lceil\frac{2^{169}}{5^{32}}\right\rceil\times\frac{1}{2^{169}}\right\rfloor=\left\lfloor\frac{x}{10^{32}}\right\rfloor,\quad 2^{94}\lt\left\lceil\frac{2^{169}}{5^{32}}\right\rceil\lt 2^{95}$

マジックナンバーがそれぞれ64bit整数を超え、 $x$ の値も64bit整数を超えているため、この計算を行うには若干の面倒が生じます。

### c. 丁度か $1$ 大きい値を求める ($d:$ 除数, $x:$ 被除数, $s,t:$ パラメータ [ $x$ の範囲に影響 ])

$\displaystyle d,s,t,x\in\mathbb{Z},\quad 0\le \min(s,t,x),\quad 2^s\le d\le 2^{(s+t)},\quad \left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)(x+2^s-1)\lt 2^{(s+t)}(d-2^s+2)$

$\displaystyle\quad\longrightarrow\quad\left\lfloor\frac{x}{d}\right\rfloor\le\left\lfloor\left\lceil\frac{x}{2^s}\right\rceil\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\frac{1}{2^t}\right\rfloor\le\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)$

- 前計算: $\displaystyle u=2^s-1,\quad m=\left\lceil\frac{2^{(s+t)}}{d}\right\rceil$
- 擬似コード: `(((x + u) >> s) * m) >> t`
- $s=0$ の場合の擬似コード: `(x * m) >> t`
- $s=0, t=W$ の場合の擬似コード: `mulh(m, x)` ※ $d=1$ の場合、 $m=\lceil 2^W/d\rceil=2^W$ となる事に注意

$s=0$ の場合のみの実用を想定。 ($s\gt 0$ だと b. の手法に比べて利点が薄いため)

$s=0, t=W=64$ の場合の利用例: ac-library の barrett reduction 実装における、除算の商の近似値を求める部分

$1\leq d\lt 2^{64}:\text{ 除数},\quad 0\leq x\lt 2^{64}:\text{ 被除数}$ $\displaystyle\quad\longrightarrow\quad\left\lfloor\frac{x}{d}\right\rfloor\leq\left\lfloor\left\lceil\frac{2^{64}}{d}\right\rceil\frac{x}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{d}\right\rfloor+1$

https://github.com/atcoder/ac-library/blob/master/atcoder/internal_math.hpp#L51-L59

#### 導出

$\displaystyle d,s,t,x\in\mathbb{Z},\quad 0\le \min(s,t,x),\quad 2^s\le d\le 2^{(s+t)},\quad \left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)(x+2^s-1)\lt 2^{(s+t)}(d-2^s+2)$

$\displaystyle\quad\longrightarrow\quad 0\leq\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)x\lt 2^{(s+t)}(d+1)-(2^s-1)\left(2^{(s+t)}+\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)\right)$

$\displaystyle\quad\longrightarrow\quad 0\leq \frac{\left(-2^{(s+t)}\right)\operatorname{mod}d}{2^{(s+t)}d}x\lt \left(\frac{d+1}{d}-\frac{(2^s-1)\left(2^{(s+t)}+\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)\right)}{2^{(s+t)}d}\right) \quad\quad$ $\displaystyle\left(\frac{1}{2^{(s+t)}d}\text{ を掛ける}\right)$

$\displaystyle\quad\longrightarrow\quad\frac{x}{d}\leq\left(\frac{x}{d}+\frac{2^s-1}{2^{(s+t)}}\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\right)\leq\left(\frac{x+2^s-1}{2^s}\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\frac{1}{2^t}\right)=\frac{\left(x+2^s-1\right)\left(2^{(s+t)}+\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)\right)}{2^{(s+t)}d}\lt\left(\frac{x}{d}+\frac{d+1}{d}\right)\leq\left(\left\lfloor\frac{x}{d}\right\rfloor+2\right)$
$\displaystyle\left(\left(\frac{x}{d}+\frac{2^s-1}{2^{(s+t)}}\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\right)=\frac{2^{(s+t)}x+(2^s-1)\left(2^{(s+t)}+\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)\right)}{2^{(s+t)}d}\text{ を足す},\ \frac{2^{(s+t)}+\left(\left(-2^{(s+t)}\right)\operatorname{mod}d\right)}{d}=\left\lceil\frac{2^{(s+t)}}{d}\right\rceil,\quad\left(\frac{x}{d}+\frac{d+1}{d}\right)\leq\left(\left\lfloor\frac{x}{d}\right\rfloor+2\right)\right)$

$\displaystyle\quad\longrightarrow\quad\frac{x}{d}\leq\left(\left\lceil\frac{x}{2^s}\right\rceil\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\frac{1}{2^t}\right)\leq\left(\frac{x+2^s-1}{2^s}\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\frac{1}{2^t}\right)\lt\left(\left\lfloor\frac{x}{d}\right\rfloor+2\right)\quad\quad$ $\displaystyle\left(\frac{x}{2^s}\leq\left\lceil\frac{x}{2^s}\right\rceil\leq\frac{x+2^s-1}{2^s}\right)$

$\displaystyle\quad\longrightarrow\quad\left\lfloor\frac{x}{d}\right\rfloor\le\left\lfloor\left\lceil\frac{x}{2^s}\right\rceil\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\frac{1}{2^t}\right\rfloor\le\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)\quad\quad$ $\left(\text{床関数(floor)を適用}\right)$

#### 例

$\displaystyle x\in\mathbb{Z},\quad 0\leq x\lt 2^{67}\lt\frac{2^{64}\times 998244353}{(-2^{64})\operatorname{mod}998244353}\quad\longrightarrow\quad\left\lfloor\frac{x}{998244353}\right\rfloor\leq\left\lfloor\left\lceil\frac{2^{64}}{998244353}\right\rceil\times\frac{x}{2^{64}}\right\rfloor=\left\lfloor\frac{2^{64}+66192443}{998244353}\times\frac{x}{2^{64}}\right\rfloor\leq\left\lfloor\frac{x}{998244353}\right\rfloor+1$

## コード例1: 除数が定数だった場合の除算/剰余算のコンパイラ最適化と等価コードの例

### コード例1a: C++

`(a * b) >> 64`, `x / 998244353`, `x mod 998244353`, `x / 7`, `x mod 7` を計算する関数 の例 (C++): https://godbolt.org/z/6fbWGsEhr
```cpp
#include <cstdint>
#ifdef _MSC_VER
#include <intrin.h>
#endif
// u64mulh_shim(a, b) == floor(a * b / 2**64)
uint64_t u64mulh_shim(uint64_t a, uint64_t b) {
#ifdef _MSC_VER // msvc https://learn.microsoft.com/ja-jp/cpp/intrinsics/umulh?view=msvc-170
    return __umulh(a, b); // Return the high 64 bits of the product of two 64-bit unsigned integers
#else // gcc/clang (128bit Integers) https://gcc.gnu.org/onlinedocs/gcc/_005f_005fint128.html
    return (uint64_t)((unsigned __int128)a * (unsigned __int128)b >> 64);
#endif
}
// div998244353_shim(x) == floor(x / 998244353)
uint64_t div998244353_shim(uint64_t x) {
    // 0 <= x < (2**93 / ((-2**93) mod 998244353)) --> floor(ceil(2**93 / 998244353) * x / 2**93) == floor(x / 998244353)
    // ceil(2**93 / 998244353) == floor((2**93 - 1) / 998244353) + 1 == 9920937979283557439
    // floor(x / 998244353) == floor((x * 9920937979283557439) / 2**93) == floor(floor((x * 9920937979283557439) / 2**64) / 2**29)
    return u64mulh_shim(x, 9920937979283557439ULL) >> 29;
}
uint64_t div998244353(uint64_t x) {
    return x / 998244353ULL;
}
// rem998244353_shim(x) == (x mod 998244353)
uint64_t rem998244353_shim(uint64_t x) {
    return x - div998244353_shim(x) * 998244353ULL;
}
uint64_t rem998244353(uint64_t x) {
    return x % 998244353ULL;
}
// div7_shim(x) == floor(x / 7)
uint64_t div7_shim(uint64_t x) {
    // 0 <= x < (2**67 / ((-2**67) mod 7)) --> floor(ceil(2**67 / 7) * x / 2**67) == floor(x / 7)
    // y := floor((ceil(2**67 / 7) - 2**64) * x / 2**64)
    // floor(x / 7) == floor((x + y) / 2**3) == floor((floor((x - y) / 2) + y) / 2**2)
    // ceil(2**67 / 7) - 2**64 == 2635249153387078803
    uint64_t y = u64mulh_shim(x, 2635249153387078803ULL);
    return (((x - y) >> 1) + y) >> 2;
}
uint64_t div7(uint64_t x) {
    return x / 7;
}
// rem7_shim(x) == (x mod 7)
uint64_t rem7_shim(uint64_t x) {
    return x - div7_shim(x) * 7ULL;
}
uint64_t rem7(uint64_t x) {
    return x % 7ULL;
}
```

#### コンパイル結果の例 (x86-64 gcc)

```nasm
# x86-64 gcc
u64mulh_shim(unsigned long, unsigned long):
    mov     rax, rdi
    mul     rsi
    mov     rax, rdx
    ret
div998244353_shim(unsigned long):
    movabs  rax, -8525806094425994177
    mul     rdi
    mov     rax, rdx
    shr     rax, 29
    ret
div998244353(unsigned long):
    movabs  rax, -8525806094425994177
    mul     rdi
    mov     rax, rdx
    shr     rax, 29
    ret
rem998244353_shim(unsigned long):
    movabs  rax, -8525806094425994177
    mul     rdi
    mov     rax, rdi
    shr     rdx, 29
    imul    rdx, rdx, 998244353
    sub     rax, rdx
    ret
rem998244353(unsigned long):
    movabs  rax, -8525806094425994177
    mul     rdi
    mov     rax, rdx
    shr     rax, 29
    imul    rdx, rax, 998244353
    mov     rax, rdi
    sub     rax, rdx
    ret
div7_shim(unsigned long):
    movabs  rax, 2635249153387078803
    mul     rdi
    sub     rdi, rdx
    shr     rdi
    lea     rax, [rdi+rdx]
    shr     rax, 2
    ret
div7(unsigned long):
    movabs  rax, 2635249153387078803
    mul     rdi
    sub     rdi, rdx
    shr     rdi
    lea     rax, [rdx+rdi]
    shr     rax, 2
    ret
rem7_shim(unsigned long):
    movabs  rax, 2635249153387078803
    mov     rcx, rdi
    mul     rdi
    mov     rax, rdi
    sub     rcx, rdx
    shr     rcx
    add     rdx, rcx
    shr     rdx, 2
    lea     rcx, [0+rdx*8]
    sub     rcx, rdx
    sub     rax, rcx
    ret
rem7(unsigned long):
    movabs  rax, 2635249153387078803
    mul     rdi
    mov     rax, rdi
    sub     rax, rdx
    shr     rax
    add     rax, rdx
    shr     rax, 2
    lea     rdx, [0+rax*8]
    sub     rdx, rax
    mov     rax, rdi
    sub     rax, rdx
    ret
```

#### コンパイル結果の例 (ARM64 gcc)

```nasm
# ARM64 gcc
u64mulh_shim(unsigned long, unsigned long):
    umulh   x0, x0, x1
    ret
div998244353_shim(unsigned long):
    mov     x1, 52287
    movk    x1, 0x5de0, lsl 16
    movk    x1, 0x4087, lsl 32
    movk    x1, 0x89ae, lsl 48
    umulh   x0, x0, x1
    lsr     x0, x0, 29
    ret
div998244353(unsigned long):
    mov     x1, 52287
    movk    x1, 0x5de0, lsl 16
    movk    x1, 0x4087, lsl 32
    movk    x1, 0x89ae, lsl 48
    umulh   x0, x0, x1
    lsr     x0, x0, 29
    ret
rem998244353_shim(unsigned long):
    mov     x2, 52287
    movk    x2, 0x5de0, lsl 16
    movk    x2, 0x4087, lsl 32
    movk    x2, 0x89ae, lsl 48
    umulh   x2, x0, x2
    lsr     x2, x2, 29
    lsl     x1, x2, 4
    sub     x1, x1, x2
    lsl     x1, x1, 3
    sub     x1, x1, x2
    add     x1, x2, x1, lsl 23
    sub     x0, x0, x1
    ret
rem998244353(unsigned long):
    mov     x2, 52287
    movk    x2, 0x5de0, lsl 16
    movk    x2, 0x4087, lsl 32
    movk    x2, 0x89ae, lsl 48
    umulh   x2, x0, x2
    lsr     x2, x2, 29
    lsl     x1, x2, 4
    sub     x1, x1, x2
    lsl     x1, x1, 3
    sub     x1, x1, x2
    add     x1, x2, x1, lsl 23
    sub     x0, x0, x1
    ret
div7_shim(unsigned long):
    mov     x1, 9363
    movk    x1, 0x9249, lsl 16
    movk    x1, 0x4924, lsl 32
    movk    x1, 0x2492, lsl 48
    umulh   x1, x0, x1
    sub     x0, x0, x1
    add     x0, x1, x0, lsr 1
    lsr     x0, x0, 2
    ret
div7(unsigned long):
    mov     x1, 9363
    movk    x1, 0x9249, lsl 16
    movk    x1, 0x4924, lsl 32
    movk    x1, 0x2492, lsl 48
    umulh   x1, x0, x1
    sub     x0, x0, x1
    add     x0, x1, x0, lsr 1
    lsr     x0, x0, 2
    ret
rem7_shim(unsigned long):
    mov     x2, 9363
    movk    x2, 0x9249, lsl 16
    movk    x2, 0x4924, lsl 32
    movk    x2, 0x2492, lsl 48
    umulh   x2, x0, x2
    sub     x1, x0, x2
    add     x1, x2, x1, lsr 1
    lsr     x1, x1, 2
    lsl     x2, x1, 3
    sub     x1, x2, x1
    sub     x0, x0, x1
    ret
rem7(unsigned long):
    mov     x2, 9363
    movk    x2, 0x9249, lsl 16
    movk    x2, 0x4924, lsl 32
    movk    x2, 0x2492, lsl 48
    umulh   x2, x0, x2
    sub     x1, x0, x2
    add     x1, x2, x1, lsr 1
    lsr     x1, x1, 2
    lsl     x2, x1, 3
    sub     x1, x2, x1
    sub     x0, x0, x1
    ret
```

### コード例1b: Rust

`(a * b) >> 64`, `x / 998244353`, `x mod 998244353`, `x / 7`, `x mod 7` を計算する関数の例 (Rust):
https://rust.godbolt.org/z/45GWzoT7K

```rust
// u64mulh_shim(a, b) == floor(a * b / 2**64)
pub fn u64mulh_shim(a: u64, b: u64) -> u64 {
    ((a as u128) * (b as u128) >> 64) as u64
}
// div998244353_shim(x) == floor(x / 998244353)
pub fn div998244353_shim(x: u64) -> u64 {
    // 0 <= x < (2**93 / ((-2**93) mod 998244353)) --> floor(ceil(2**93 / 998244353) * x / 2**93) == floor(x / 998244353)
    // ceil(2**93 / 998244353) == floor((2**93 - 1) / 998244353) + 1 == 9920937979283557439
    // floor(x / 998244353) == floor((x * 9920937979283557439) / 2**93) == floor(floor((x * 9920937979283557439) / 2**64) / 2**29)
    u64mulh_shim(x, 9920937979283557439) >> 29
}
pub fn div998244353(x: u64) -> u64 {
    x / 998244353
}
// rem998244353_shim(x) == (x mod 998244353)
pub fn rem998244353_shim(x: u64) -> u64 {
    x - div998244353_shim(x) * 998244353
}
pub fn rem998244353(x: u64) -> u64 {
    x % 998244353
}
// div7_shim(x) == floor(x / 7)
pub fn div7_shim(x: u64) -> u64 {
    // 0 <= x < (2**67 / ((-2**67) mod 7)) --> floor(ceil(2**67 / 7) * x / 2**67) == floor(x / 7)
    // y := floor((ceil(2**67 / 7) - 2**64) * x / 2**64)
    // floor(x / 7) == floor((x + y) / 2**3) == floor((floor((x - y) / 2) + y) / 2**2)
    // ceil(2**67 / 7) - 2**64 == 2635249153387078803
    let y = u64mulh_shim(x, 2635249153387078803);
    (((x - y) >> 1) + y) >> 2
}
pub fn div7(x: u64) -> u64 {
    x / 7
}
// rem7_shim(x) == (x mod 7)
pub fn rem7_shim(x: u64) -> u64 {
    x - div7_shim(x) * 7
}
pub fn rem7(x: u64) -> u64 {
    x % 7
}
```

#### コンパイル結果の例 (rustc x86_64-unknown-linux-gnu)

```nasm
# rustc x86_64-unknown-linux-gnu
u64mulh_shim:
    mov     rax, rsi
    mul     rdi
    mov     rax, rdx
    ret

div998244353_shim:
    mov     rax, rdi
    movabs  rcx, -8525806094425994177
    mul     rcx
    mov     rax, rdx
    shr     rax, 29
    ret

div998244353:
    mov     rax, rdi
    movabs  rcx, -8525806094425994177
    mul     rcx
    mov     rax, rdx
    shr     rax, 29
    ret

rem998244353_shim:
    movabs  rcx, -8525806094425994177
    mov     rax, rdi
    mul     rcx
    shr     rdx, 29
    imul    rax, rdx, -998244353
    add     rax, rdi
    ret

rem998244353:
    movabs  rcx, -8525806094425994177
    mov     rax, rdi
    mul     rcx
    shr     rdx, 29
    imul    rax, rdx, 998244353
    sub     rdi, rax
    mov     rax, rdi
    ret

div7_shim:
    movabs  rcx, 2635249153387078803
    mov     rax, rdi
    mul     rcx
    sub     rdi, rdx
    shr     rdi
    lea     rax, [rdi + rdx]
    shr     rax, 2
    ret

div7:
    movabs  rcx, 2635249153387078803
    mov     rax, rdi
    mul     rcx
    sub     rdi, rdx
    shr     rdi
    lea     rax, [rdi + rdx]
    shr     rax, 2
    ret

rem7_shim:
    movabs  rcx, 2635249153387078803
    mov     rax, rdi
    mul     rcx
    mov     rax, rdi
    sub     rax, rdx
    shr     rax
    add     rax, rdx
    shr     rax, 2
    lea     rcx, [8*rax]
    sub     rax, rcx
    add     rax, rdi
    ret

rem7:
    movabs  rcx, 2635249153387078803
    mov     rax, rdi
    mul     rcx
    mov     rax, rdi
    sub     rax, rdx
    shr     rax
    add     rax, rdx
    shr     rax, 2
    lea     rcx, [8*rax]
    sub     rax, rcx
    add     rax, rdi
    ret
```

#### コンパイル結果の例 (rustc aarch64-unknown-linux-gnu)

```nasm
# rustc aarch64-unknown-linux-gnu
u64mulh_shim:
    umulh   x0, x1, x0
    ret

div998244353_shim:
    mov     x8, #52287
    movk    x8, #24032, lsl #16
    movk    x8, #16519, lsl #32
    movk    x8, #35246, lsl #48
    umulh   x8, x0, x8
    lsr     x0, x8, #29
    ret

div998244353:
    mov     x8, #52287
    movk    x8, #24032, lsl #16
    movk    x8, #16519, lsl #32
    movk    x8, #35246, lsl #48
    umulh   x8, x0, x8
    lsr     x0, x8, #29
    ret

rem998244353_shim:
    mov     x8, #52287
    mov     x9, #-998244353
    movk    x8, #24032, lsl #16
    movk    x8, #16519, lsl #32
    movk    x8, #35246, lsl #48
    umulh   x8, x0, x8
    lsr     x8, x8, #29
    madd    x0, x8, x9, x0
    ret

rem998244353:
    mov     x8, #52287
    mov     w9, #1
    movk    x8, #24032, lsl #16
    movk    w9, #15232, lsl #16
    movk    x8, #16519, lsl #32
    movk    x8, #35246, lsl #48
    umulh   x8, x0, x8
    lsr     x8, x8, #29
    msub    x0, x8, x9, x0
    ret

div7_shim:
    mov     x8, #9363
    movk    x8, #37449, lsl #16
    movk    x8, #18724, lsl #32
    movk    x8, #9362, lsl #48
    umulh   x8, x0, x8
    sub     x9, x0, x8
    add     x8, x8, x9, lsr #1
    lsr     x0, x8, #2
    ret

div7:
    mov     x8, #9363
    movk    x8, #37449, lsl #16
    movk    x8, #18724, lsl #32
    movk    x8, #9362, lsl #48
    umulh   x8, x0, x8
    sub     x9, x0, x8
    add     x8, x8, x9, lsr #1
    lsr     x0, x8, #2
    ret

rem7_shim:
    mov     x8, #9363
    movk    x8, #37449, lsl #16
    movk    x8, #18724, lsl #32
    movk    x8, #9362, lsl #48
    umulh   x8, x0, x8
    sub     x9, x0, x8
    add     x8, x8, x9, lsr #1
    lsr     x8, x8, #2
    sub     x8, x8, x8, lsl #3
    add     x0, x8, x0
    ret

rem7:
    mov     x8, #9363
    movk    x8, #37449, lsl #16
    movk    x8, #18724, lsl #32
    movk    x8, #9362, lsl #48
    umulh   x8, x0, x8
    sub     x9, x0, x8
    add     x8, x8, x9, lsr #1
    lsr     x8, x8, #2
    sub     x8, x8, x8, lsl #3
    add     x0, x0, x8
    ret
```

## $\mathrm{mod}(2^{61}-1)$ の計算 (rolling hash)

- C++: https://godbolt.org/z/E7xoGqEha
- Rust: https://rust.godbolt.org/z/1WGK86rda

https://yosupo.hatenablog.com/entry/2023/08/06/181942

## 導出2

- $\lfloor x\rfloor:$ 床関数 (小数点以下 $-\infty$ 方向に切り下げ)
- $\lceil x\rceil:$ 天井関数 (小数点以下 $+\infty$ 方向に切り上げ)

$\displaystyle n,d,s,t,x\in\mathbb{Z},\quad 0\leq\min(s,t,x),\quad 2^s\leq d\leq 2^{(s+t)},$

$\displaystyle n\ge 2,\quad 2^{(s+t)}/d\leq 2^{(n+1)}-1,\quad x/2^s\leq 2^n-1$

$n:$ 整数型のbit数 $,\ d:$ 除数 $,\ x:$ 被除数 $,\ s,t:$ 演算が正しく動く区間を調整するためのパラメータ

問題: 

1. $s=0$、正整数 $d$、非負整数 $t$ が与えられた時、 $\displaystyle\left\lfloor\left\lceil\frac{2^t}{d}\right\rceil\times\frac{x}{2^t}\right\rfloor=\left\lfloor\frac{x}{d}\right\rfloor$ が成り立つ 非負整数 $x$ の区間を求めよ。
2. 正整数 $d$、非負整数 $s,t$ が与えられた時、 $\displaystyle\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)\leq\left\lfloor\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\times\left\lfloor\frac{x}{2^s}\right\rfloor\times\frac{1}{2^t}\right\rfloor\leq\left\lfloor\frac{x}{d}\right\rfloor$ が成り立つ 非負整数 $x$ の区間を求めよ。
3. 正整数 $d$、非負整数 $s,t$ が与えられた時、 $\displaystyle\left\lfloor\frac{x}{d}\right\rfloor\leq\left\lfloor\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\times\left\lceil\frac{x}{2^s}\right\rceil\times\frac{1}{2^t}\right\rfloor\leq\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)$ が成り立つ 非負整数 $x$ の区間を求めよ。
4. 正整数 $d$、非負整数 $s,t$ が与えられた時、 $\displaystyle\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)\leq\left\lfloor\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\times\left\lfloor\frac{x}{2^s}\right\rfloor\times\frac{1}{2^t}\right\rfloor\leq\left\lfloor\frac{x}{d}\right\rfloor$ が成り立つ 非負整数 $x$ の区間を求めよ。

暫定解: (※) の所で極値を代入していますが、常にそうであるとは限らないため、個別のケースでは問に示した不等式が成り立つ区間が広がる事があります

1. 非負整数 $r,w$ を用いて $\displaystyle\left\lceil\frac{2^t}{d}\right\rceil=\frac{2^t+r}{d},\quad \left\lfloor\frac{x}{d}\right\rfloor=\frac{x-w}{d}$ とおくと
    > $\displaystyle 0\leq r=\left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)\lt d,\quad 0\leq w=(x\,\mathrm{mod}\,d)\lt d$
    >
    > $\displaystyle \left(\left\lceil\frac{2^t}{d}\right\rceil\times\frac{x}{2^t}\right)=\left(\frac{2^t+r}{d}\times\frac{x}{2^t}\right)=\left(\frac{d}{x}+\frac{rx}{2^td}\right)\lt\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)=\left(\frac{x}{d}+1-\frac{w}{d}\right)$
    >
    > $\displaystyle\quad\longrightarrow\quad \frac{rx}{2^td}\lt 1-\frac{w}{d}$
    >
    > $\displaystyle\quad\longrightarrow\quad rx\lt 2^t(d-w)$
    >
    > (※) 最悪ケースの $w=(d-1)$ をとると $0\leq x$ より
    >
    > $\displaystyle 0\leq ((-2^t)\,\mathrm{mod}\,d)x\lt 2^t$
2. 非負整数 $r,v,w$ を用いて $\displaystyle\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor=\frac{2^{(s+t)}-r}{d},\quad\left\lfloor\frac{x}{2^s}\right\rfloor=\frac{x-v}{2^s},\quad\left\lfloor\frac{x}{d}\right\rfloor=\frac{x-w}{d}$ とおくと
    > $\displaystyle 0\leq r=(2^{(s+t)}\,\mathrm{mod}\,d)\lt d,\quad 0\leq v=(x\,\mathrm{mod}\,2^s)\lt 2^s,\quad 0\leq w=(x\,\mathrm{mod}\,d)\lt d$
    >
    > $\displaystyle\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)=\frac{x-d-w}{d}\leq\left(\left\lfloor\frac{2^{(s+t)}}{d}\right\rfloor\times\left\lfloor\frac{x}{2^s}\right\rfloor\times\frac{1}{2^t}\right)=\frac{(2^{(s+t)}-r)(x-v)}{2^{(s+t)}d}\lt\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)=\frac{x+d-w}{d}$
    >
    > $\displaystyle\quad\longrightarrow\quad 2^{(s+t)}(x-d-w)\leq(2^{(s+t)}-r)(x-v)\lt 2^{(s+t)}(x-w+d)$
    >
    > $\displaystyle\quad\longrightarrow\quad rx\leq 2^{(s+t)}(d+w)-(2^{(s+t)}-r)v$
    >
    > (※) $2^{s+t}\geq d\gt r$ より、最悪ケースの $v=(2^s-1), w=0$ をとると
    >
    > $\displaystyle (2^{(s+t)}\,\mathrm{mod}\,d)(x-2^s+1)\leq 2^{(s+t)}(d-2^s+1)$
3. 非負整数 $r,v,w$ を用いて $\displaystyle\left\lceil\frac{2^{(s+t)}}{d}\right\rceil=\frac{2^{(s+t)}+r}{d},\quad\left\lceil\frac{x}{2^s}\right\rceil=\frac{x+v}{2^s},\quad\left\lfloor\frac{x}{d}\right\rfloor=\frac{x-w}{d}$ とおくと
    > $\displaystyle 0\leq r=\left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)\lt d,\quad 0\leq v=((-x)\,\mathrm{mod}\,2^s)\lt 2^s,\quad 0\leq w=(x\,\mathrm{mod}\,d)\lt d$
    >
    > $\displaystyle\left\lfloor\frac{x}{d}\right\rfloor=\frac{x-w}{d}\leq\left(\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\times\left\lceil\frac{x}{2^s}\right\rceil\times\frac{1}{2^t}\right)=\frac{(2^{(s+t)}+r)(x+v)}{2^{(s+t)}d}\lt\left(\left\lfloor\frac{x}{d}\right\rfloor+2\right)=\frac{x+2d-w}{d}$
    >
    > $\displaystyle\quad\longrightarrow\quad 0\leq\frac{(2^{(s+t)}+r)(x+v)-2^{(s+t)}(x-w)}{2^{(s+t)}d}\lt 2$
    >
    > $\displaystyle\quad\longrightarrow\quad rx\lt 2^{(s+t)}\times 2d-\left((2^{(s+t)}+r)v+2^{(s+t)}w\right)$
    >
    > (※) 最悪ケースの $v=(2^s-1), w=(d-1)$ をとると
    >
    > $\displaystyle\left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)(x+2^s-1)\lt 2^{(s+t)}(d-2^s+2)$
4. 非負整数 $r,v,w$ を用いて $\displaystyle\left\lceil\frac{2^{(s+t)}}{d}\right\rceil=\frac{2^{(s+t)}+r}{d},\quad\left\lfloor\frac{x}{2^s}\right\rfloor=\frac{x-v}{2^s},\quad\left\lfloor\frac{x}{d}\right\rfloor=\frac{x-w}{d}$ とおくと
    > $\displaystyle 0\leq r=\left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)\lt d,\quad 0\leq v=(x\,\mathrm{mod}\,2^s)\lt 2^s,\quad 0\leq w=(x\,\mathrm{mod}\,d)\lt d$
    >
    > $\displaystyle\left(\left\lfloor\frac{x}{d}\right\rfloor-1\right)=\left(\frac{x-w}{d}-1\right)\leq\left(\left\lceil\frac{2^{(s+t)}}{d}\right\rceil\times\left\lfloor\frac{x}{2^s}\right\rfloor\times\frac{1}{2^t}\right)=\frac{(2^{(s+t)}+r)(x-v)}{2^{(s+t)}d}\lt\left(\left\lfloor\frac{x}{d}\right\rfloor+1\right)=\left(\frac{x-w}{d}+1\right)$
    >
    > $\displaystyle\quad\longrightarrow\quad -1\leq\frac{(2^{s+t}+r)(x-v)-2^{(s+t)}(x-w)}{2^{(s+t)}d}\lt 1$
    >
    > $\displaystyle\quad\longrightarrow\quad -2^{(s+t)}d\leq rx-(2^{(s+t)}+r)v+2^{(s+t)}w\lt 2^{(s+t)}d$
    >
    > (※) 左側の不等式には最悪ケースの $v=2^s-1,w=0$ を、 右側の不等式には最悪ケースの $v=0,w=d-1$ をとると
    >
    > $\displaystyle \left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)(2^s-1)-2^{(s+t)}(d-2^s+1)\leq \left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)x\lt 2^{(s+t)}$
    
$\displaystyle \left(\left(-2^{(s+t)}\right)\,\mathrm{mod}\,d\right)(2^s-1)\leq 2^{(s+t)}(d-2^s+1)$
