---
title: "モジュラ逆数の計算法"
emoji: "🦀"
type: "tech"
topics:
  - "rust"
  - "数学"
  - "競技プログラミング"
  - "アルゴリズム"
published: true
published_at: "2023-04-20 17:02"
---

# モジュラ逆数の計算

https://ja.wikipedia.org/wiki/%E3%83%A2%E3%82%B8%E3%83%A5%E3%83%A9%E9%80%86%E6%95%B0

$0\lt n \lt m\le 2^{64}$ 辺りの範囲で $nn^{-1}\equiv 1\ (\operatorname{mod}m)$ が成り立つ モジュラ逆数 $n^{-1}\mod m$ を求める方法を考えてみます。

(サンプルコードでは主に $0\lt n\lt m\lt 2^{32}$ を前提として書きます)

モジュラ逆数 $n^{-1}\mod m$ が存在するための必要十分条件は、 $n$ と $m$ の 最大公約数 $\gcd(n,m)$ が $\gcd(n,m)=1$ を満たすことです。

ここでは、3つの方法でのモジュラ逆数の求め方を紹介します。

- オイラーの定理を用いる方法 ( $m$ が素数である場合、 $m$ の素因数分解が既知である場合 )
- 変則的な拡張ユークリッド互除法を用いる方法 ( 一般の自然数 $m$ に対する場合 )
- ニュートン・ラフソン法を用いる方法 ( $m$ が累乗数である場合 )

# $m$ が素数、もしくは $m$ の素因数分解が既知である場合のモジュラ逆数 (オイラーの定理)

オイラーのトーシェント関数 $\varphi(m)$ を用いると、オイラーの定理より

- $\gcd(n,m)=1\quad\Longrightarrow\quad n^{\varphi(m)}\equiv 1\ (\operatorname{mod}m)$
- $\gcd(n,m)=1\quad\Longrightarrow\quad n^{-1}\equiv n^{\varphi(m)-1}\ (\operatorname{mod}m)$

$n^{\varphi(m)-1}\ (\operatorname{mod}m)$ は、繰り返し二乗法（バイナリ法）などを使って計算できます。

$m$ が素数であれば、$\varphi(m)=m-1$ であるので、 $\gcd(n,m)=1$ ならば $n^{-1}\equiv n^{m-2}\ (\operatorname{mod}m)$ と計算できます。

$m$ が素数でなければ、 $\varphi(m)$ の値は以下のように計算する必要があります。

## オイラーのトーシェント関数

オイラーのトーシェント関数 $\varphi(m)$ は、正の整数 $m$ に対して $m$ と互いに素である $1$ 以上 $m$ 以下の自然数の個数です。

$$\displaystyle\varphi(m)=\sum_{1\le n\le m\atop \gcd(n,m)}1$$

- $m$ が素数であれば、 $\varphi(m)=m-1$
- $m$ が素数ではなく、 $m$ の素因数分解が以下のように ( $p_1,p_2,p_3,\dots,p_d$ は互いに異なる素数 ) 既知であれば、

$$\displaystyle m=\prod_{k=1}^d p_k^{\,e_k}=p_1^{\,e_1}\times p_2^{\,e_2}\times p_3^{\,e_3}\times\cdots\times p_d^{\,e_d}$$


$$
\varphi(m)=m\prod_{k=1}^d\left(1-\frac{1}{p_k}\right)=m\left(1-\frac{1}{p_1}\right)\left(1-\frac{1}{p_2}\right)\left(1-\frac{1}{p_3}\right)\cdots\left(1-\frac{1}{p_d}\right)
$$

## サンプルコード (Rust)

```rust
// モジュラ逆数 $n^{-1} \mod 998244353$ をオイラーの定理から計算 (998244353 は素数)
pub fn invmod_euler_998244353(n: u32) -> Option<u32> {
    const MOD: u32 = 998_244_353;
    let nm = n % MOD;
    if nm == 0 {
        // n が m の倍数なら gcd(n,m) != 1 なので、モジュラ逆数は存在しない
        return None;
    }
    // 内部計算には64bit整数を使用(変数r,変数s)
    let (mut r, mut s, mut t) = (1u64, u64::from(nm), MOD - 2);
    // 繰り返し二乗法（バイナリ法）にて $n^{M-2}$ を計算
    while t > 0 {
        if t % 2 == 1 {
            r = (r * s) % u64::from(MOD);
        }
        s = (s * s) % u64::from(MOD);
        t /= 2;
    }
    // 32ビット整数に変換して出力
    Some(r as u32)
}
```

# $m$ が累乗数の場合のモジュラ逆数 (ニュートン・ラフソン法)

自然数 $p,q$ を用いて $m=p^w$ とできる場合、自然数 $r$ と ニュートン・ラフソン法における、実数の逆数近似計算の漸化式 $x_{i+1}=x_i(2-nx_i)$ を用いて

$nx_i\equiv 1\ (\operatorname{mod}p^r),\quad x_{i+1}=x_i(2-nx_i) \quad\Longrightarrow\quad nx_{i+1}\equiv 1\ (\operatorname{mod}p^{2r})$

と計算する事ができます。また、 $w\le 2r$ であれば、

$w\le 2r$ かつ $nx_i\equiv 1\ (\operatorname{mod}p^r) \quad\Longrightarrow\quad nx_{i+1}\equiv 1\ (\operatorname{mod}p^w)$

となります。

証明:

$nx_i=1+kp^r$ とおくと、

$nx_{i+1}=nx_i(2-nx_i)=(1+kp^r)(2-(1+kp^r))$

$=(1+kp^r)(1-kp^r)=(1-k^2p^{2r})\equiv 1\ (\operatorname{mod}p^{2r})$

また、 1回の反復で有効桁数を2倍ではなく3倍にする方法もあります。

$nx_i\equiv 1\ (\operatorname{mod}p^r),\quad x_{i+1}=x_i(nx_i(nx_i-3)+3)\quad\Longrightarrow\quad nx_{i+1}\equiv 1\ (\operatorname{mod}p^{3r})$

証明:

$nx_i=1+kp^r$ とおくと、
$nx_{i+1}=nx_i(nx_i(nx_i-3)+3)$
$=1+(nx_i-1)^3=1+k^3p^{3r}\equiv 1\quad(\operatorname{mod}\ p^{3r})$

## 例: $m$ が $2$ の累乗数 ( $m=2^w$ ) の場合

$n$ が奇数 : $n\equiv 1\ (\operatorname{mod}2)$ であれば、 $n^2\equiv 1\ (\operatorname{mod}8)$ となります。

> - 奇数 $n$ を 整数 $k$ を用いて $n=2k+1$ とすると、 $\displaystyle n^2=8\times\left(\frac{k(k+1)}{2}\right)+1$
> - $k(k+1)$ は偶数なので、 $\displaystyle \left(\frac{k(k+1)}{2}\right)$ は整数
> - よって、 $n^2=8\times(\text{整数})+1$ より $n^2\equiv 1\ (\operatorname{mod}8)$

また、 $n$ が奇数 : $n\equiv 1\ (\operatorname{mod}2)$ であれば、 $(n\times \operatorname{XOR}(3n,2))\equiv 1\ (\operatorname{mod}32)$ となります。

> - $\operatorname{XOR}$ は排他的論理和によるビット演算関数。
> - $\displaystyle\operatorname{XOR}(3n,2)=\begin{cases}3n-2&(n\equiv 1\ (\operatorname{mod}4)\ \Leftrightarrow\ 3n\equiv 3\ (\operatorname{mod}4))\\3n+2&(n\equiv 3\ (\operatorname{mod}4)\ \Leftrightarrow\ 3n\equiv 1\ (\operatorname{mod}4))\end{cases}$
> - $n\equiv 1\ (\operatorname{mod}4)$ のとき、 $n=4k+1$ とすると 
>   $$\begin{align*}&\quad n\times\operatorname{XOR}(3n, 2)\\&=(4k+1)(3(4k+1)-2)\\ &=32\left(k^2+\frac{k(k+1)}{2}\right)+1\\ &\equiv 1\ (\operatorname{mod}32)\end{align*}$$
>   - $k(k+1)$ は偶数なので、 $\displaystyle\left(k^2+\frac{k(k+1)}{2}\right)$ は整数
> - $n\equiv 3\ (\operatorname{mod}4)$ のとき、 $n=4k-1$ とすると 
>   $$\begin{align*}&\quad n\times\operatorname{XOR}(3n, 2)\\&=(4k-1)(3(4k-1)+2)\\ &=32\left(k^2+\frac{k(k-1)}{2}\right)+1\\ &\equiv 1\ (\operatorname{mod}32)\end{align*}$$
>   - $k(k-1)$ は偶数なので、 $\displaystyle\left(k^2+\frac{k(k-1)}{2}\right)$ は整数
> - よって、 $n$ が奇数 : $n\equiv 1\ (\operatorname{mod}2)$ であれば、 $(n\times \operatorname{XOR}(3n,2))\equiv 1\ (\operatorname{mod}32)$

$m=2^w$ とすると、

$x_0=n$ を初期値とした場合、

1. $x_0=n\quad\Longrightarrow\quad nx_0\equiv 1\ (\operatorname{mod}2^{\min(3,w)})$
2. $x_1=x_0(2-nx_0)\mod 2^w\quad\Longrightarrow\quad nx_1\equiv 1\ (\operatorname{mod}2^{\min(6,w)})$
3. $x_2=x_1(2-nx_1)\mod 2^w\quad\Longrightarrow\quad nx_2\equiv 1\ (\operatorname{mod}2^{\min(12,w)})$
4. $x_3=x_2(2-nx_2)\mod 2^w\quad\Longrightarrow\quad nx_3\equiv 1\ (\operatorname{mod}2^{\min(24,w)})$
5. $x_4=x_3(2-nx_3)\mod 2^w\quad\Longrightarrow\quad nx_4\equiv 1\ (\operatorname{mod}2^{\min(48,w)})$
5. $x_5=x_4(2-nx_4)\mod 2^w\quad\Longrightarrow\quad nx_5\equiv 1\ (\operatorname{mod}2^{\min(96,w)})$

$x_0=\operatorname{XOR}(3n,2)\mod 2^w$ を初期値とした場合、

1. $x_0=\operatorname{XOR}(3n,2)\mod 2^w\quad\Longrightarrow\quad nx_0\equiv 1\ (\operatorname{mod}2^{\min(5,w)})$
2. $x_1=x_0(2-nx_0)\mod 2^w\quad\Longrightarrow\quad nx_1\equiv 1\ (\operatorname{mod}2^{\min(10,w)})$
3. $x_2=x_1(2-nx_1)\mod 2^w\quad\Longrightarrow\quad nx_2\equiv 1\ (\operatorname{mod}2^{\min(20,w)})$
4. $x_3=x_2(2-nx_2)\mod 2^w\quad\Longrightarrow\quad nx_3\equiv 1\ (\operatorname{mod}2^{\min(40,w)})$
5. $x_4=x_3(2-nx_3)\mod 2^w\quad\Longrightarrow\quad nx_4\equiv 1\ (\operatorname{mod}2^{\min(80,w)})$

のように $x_{i+1}=x_i(2-nx_i)\mod 2^w$ の漸化式を必要な回数 ( 前者 $x_0=n$ なら $\max(\lceil\log_2(w / 3)\rceil,0)$ 回、 後者 $x_0=\operatorname{XOR}(3n,2)\mod 2^w$ なら $\max(\lceil\log_2(w / 5)\rceil,0)$ 回 ) 計算すれば、 $n^{-1}\equiv x_{\max(\lceil\log_2(w/k)\rceil,0)}\ (\operatorname{mod}2^w)$ を計算することができます。

$\operatorname{mod}2^{64}$ におけるモジュラ逆数 ( $w=64$ ) を求めるのに、 加減算・ビットXOR・3倍 ( x86_64 CPUの場合、 `tmp = n * 3` の計算を乗算を使わず `tmp = n + (n << 1)` のように変形して1サイクル以下で実行できる ) の計算に 1サイクル、乗算に 3サイクル の実行時間が掛かるとした実行時間の例 (30サイクル)。

```
 0. `tmp0 := n * 3 mod 2^{64}` : 1cycle
 1. `x0 := xor(tmp0, 2)` : 1cycle
 2. `tmp1 := n * x0 mod 2^{64}` : 3cycle
 3.
 4.
 5. `tmp2 := 2 - tmp1 mod 2^{64}` : 1cycle
 6. `x1 := x0 * tmp2 mod 2^{64}` : 3cycle
 7.
 8.
 9. `tmp3 := n * x1 mod 2^{64}` : 3cycle
10.
11.
12. `tmp4 := 2 - tmp3 mod 2^{64}` : 1cycle
13. `x2 := x1 * tmp4 mod 2^{64}` : 3cycle
14.
15.
16. `tmp5 := n * x2 mod 2^{64}` : 3cycle
17.
18.
19. `tmp6 := 2 - tmp5 mod 2^{64}` : 1cycle
20. `x3 := x1 * tmp6 mod 2^{64}` : 3cycle
21.
22.
23. `tmp7 := n * x3 mod 2^{64}` : 3cycle
24.
25.
26. `tmp8 := 2 - tmp7 mod 2^{64}` : 1cycle
27. `x4 := x1 * tmp8 mod 2^{64}` : 3cycle
28.
29.
30. `return x4`
```

## Jeffrey Hurchalla’s method

さらにこの漸化式 $x_{i+1}=x_i(2-nx_i)\mod 2^w$ を変形し、命令レベルの並列性（Instruction-level parallelism、ILP）を向上させて実行時間の短縮を図る手法が [An Improved Integer Modular Multiplicative Inverse (modulo $2^w$) arxiv.org/abs/2204.04342](https://arxiv.org/abs/2204.04342) として発表されています。

https://ja.wikipedia.org/wiki/%E5%91%BD%E4%BB%A4%E3%83%AC%E3%83%99%E3%83%AB%E3%81%AE%E4%B8%A6%E5%88%97%E6%80%A7

https://arxiv.org/abs/2204.04342

1. $x_0=\operatorname{XOR}(3n,2)\mod 2^w,\quad y_0=1-nx_0\mod 2^w\quad\Longrightarrow\quad nx_0\equiv 1\ (\operatorname{mod}2^{\min(5,w)})$
2. $x_1=x_0(1+y_0)\mod 2^w,\quad y_1=(y_0)^2\mod 2^w\quad\Longrightarrow\quad nx_1\equiv 1\ (\operatorname{mod}2^{\min(10,w)})$
3. $x_2=x_1(1+y_1)\mod 2^w,\quad y_2=(y_1)^2\mod 2^w\quad\Longrightarrow\quad nx_2\equiv 1\ (\operatorname{mod}2^{\min(20,w)})$
4. $x_3=x_2(1+y_2)\mod 2^w,\quad y_3=(y_2)^2\mod 2^w\quad\Longrightarrow\quad nx_3\equiv 1\ (\operatorname{mod}2^{\min(40,w)})$
5. $x_4=x_3(1+y_3)\mod 2^w,\quad y_4=(y_3)^2\mod 2^w\quad\Longrightarrow\quad nx_4\equiv 1\ (\operatorname{mod}2^{\min(80,w)})$

このように ニュートン・ラフソン法の漸化式 $x_{i+1}=x_i(2-nx_i)\mod 2^w$ を変形してみると、 $x_i$ の値については Jeffrey Hurchalla’s method との間でそれぞれ対応する関係にあるのが分かります。

- $nx_0\equiv(1-y_0)\quad(\operatorname{mod}2^w)$
- $x_1\equiv x_0(2-nx_0)\equiv x_0(1+y_0)\quad(\operatorname{mod}2^w)$
- $nx_1\equiv nx_0(1+y_0)\equiv(1-y_0)(1+y_0)\equiv(1-(y_0)^2)\quad(\operatorname{mod}2^w)$
- $x_2\equiv x_1(2-nx_1)\equiv x_1(1+(y_0)^2)\quad(\operatorname{mod}2^w)$
- $nx_2\equiv nx_1(1+(y_0)^2)\equiv(1-(y_0)^2)(1+(y_0)^2)\equiv(1-(y_0)^4)\quad(\operatorname{mod}2^w)$
- $x_3\equiv x_2(2-nx_2)\equiv x_2(1+(y_0)^4)\quad(\operatorname{mod}2^w)$
- $nx_3\equiv nx_2(1+(y_0)^4)\equiv(1-(y_0)^4)(1+(y_0)^4)\equiv(1-(y_0)^8)\quad(\operatorname{mod}2^w)$
- $x_4\equiv x_3(2-nx_3)\equiv x_3(1+(y_0)^8)\quad(\operatorname{mod}2^w)$
- $nx_4\equiv nx_3(1+(y_0)^8)\equiv(1-(y_0)^8)(1+(y_0)^8)\equiv(1-(y_0)^{16})\quad(\operatorname{mod}2^w)$

$w=64$ の時、 加減算・ビットXOR・3倍 ( x86_64 CPUの場合、 `tmp = n * 3` の計算を乗算を使わず `tmp = n + (n << 1)` のように変形して1サイクル以下で実行できる ) の計算に 1サイクル、乗算に 3サイクル の実行時間が掛かるとして、以下のように 19サイクル で $\operatorname{mod}2^{64}$ におけるモジュラ逆数を求めることができます。

```
 0. `tmp0 := n * 3 mod 2^{64}` : 1cycle
 1. `x0 := xor(tmp0, 2)` : 1cycle
 2. `tmp1 := n * x0 mod 2^{64}` : 3cycle
 3.
 4.
 5. `y0 := 1 - tmp1 mod 2^{64}` : 1cycle
 6. `tmp2 := y0 + 1 mod 2^{64}` : 1cycle, `y1 := y0 * y0 mod 2^{64}`: 3cycle
 7. `x1 := x0 * tmp2 mod 2^{64}` : 3cycle, 
 8.
 9. `tmp3 := y1 + 1 mod 2^{64}` : 1cycle, `y2 := y1 * y1 mod 2^{64}` : 3cycle
10. `x2 := x1 * tmp3 mod 2^{64}` : 3cycle
11.
12. `tmp4 := y2 + 1 mod 2^{64}` : 1cycle, `y3 := y2 * y2 mod 2^{64}` : 3cycle
13. `x3 := x2 * tmp4 mod 2^{64}` : 3cycle
14.
15. `tmp5 := y3 + 1 mod 2^{64}` : 1cycle
16. `x4 := x3 * tmp5 mod 2^{64}` : 3cycle
17.
18.
19. `return x4`
```

## サンプルコード (Rust)

```rust
// モジュラ逆数 $n^{-1} \mod 2^{32}$ の計算
pub fn invmod_u32(n: u32) -> Option<u32> {
    // n が偶数だとモジュラ逆数は存在しない
    if n % 2 == 0 {
        return None;
    }
    // Jeffrey Hurchalla’s method: https://arxiv.org/abs/2204.04342
    let x0 = n.wrapping_mul(3) ^ 2; // x0 := BITWISEXOR(n * 3, 2) mod 2^{32}
    let y0 = 1u32.wrapping_sub(n.wrapping_mul(x0)); // y0 := 1 - n * x0 mod 2^{32}
    let x1 = x0.wrapping_mul(y0.wrapping_add(1)); // x1 := x0 * (y0 + 1) mod 2^{32}
    let y1 = y0.wrapping_mul(y0); // y1 := y0 * y0 mod 2^{32}
    let x2 = x1.wrapping_mul(y1.wrapping_add(1)); // x2 := x1 * (y1 + 1) mod 2^{32}
    let y2 = y1.wrapping_mul(y1); // y2 := y1 * y1 mod 2^{32}
    let x3 = x2.wrapping_mul(y2.wrapping_add(1)); // x3 := x2 * (y2 + 1) mod 2^{32}
    Some(x3) // return x3
}
// モジュラ逆数 $n^{-1} \mod 2^{64}$ の計算
pub fn invmod_u64(n: u64) -> Option<u64> {
    // n が偶数だとモジュラ逆数は存在しない
    if n % 2 == 0 {
        return None;
    }
    // Jeffrey Hurchalla’s method: https://arxiv.org/abs/2204.04342
    let x0 = n.wrapping_mul(3) ^ 2; // x0 := BITWISEXOR(n * 3, 2) mod 2^{64}
    let y0 = 1u64.wrapping_sub(n.wrapping_mul(x0)); // y0 := 1 - n * x0 mod 2^{64}
    let x1 = x0.wrapping_mul(y0.wrapping_add(1)); // x1 := x0 * (y0 + 1) mod 2^{64}
    let y1 = y0.wrapping_mul(y0); // y1 := y0 * y0 mod 2^{64}
    let x2 = x1.wrapping_mul(y1.wrapping_add(1)); // x2 := x1 * (y1 + 1) mod 2^{64}
    let y2 = y1.wrapping_mul(y1); // y2 := y1 * y1 mod 2^{64}
    let x3 = x2.wrapping_mul(y2.wrapping_add(1)); // x3 := x2 * (y2 + 1) mod 2^{64}
    let y3 = y2.wrapping_mul(y2); // y3 := y2 * y2 mod 2^{64}
    let x4 = x3.wrapping_mul(y3.wrapping_add(1)); // x4 := x3 * (y3 + 1) mod 2^{64}
    Some(x4) // return x4
}
```

# 一般の $m$ に関するモジュラ逆数 (変則・拡張ユークリッド互除法)

拡張ユークリッドの互除法によるモジュラ逆数の計算は、多くの資料では符号付き整数型を用いた実装が多いですが、ここでは少し式を工夫して符号無し整数型を用いた実装を考えてみます。

符号付き整数型を用いた拡張ユークリッド互除法のコード例 (Rust):

```rust
// extgcd(a, b) -> (g, x, y)
// ax + by = gcd(a,b) = g を満たす整数の組 (g, x, y) をひとつ返す
pub fn extgcd(a: i64, b: i64) -> (i64, i64, i64) {
    match b {
        0 => (a, 1, 0),
        _ => {
            let (d, y, x) = extgcd(b, a % b);
            (d, x, y - (a / b) * x)
        },
    }
}
```

## 対象

このような形の連立方程式を変形する操作を考えてみます。

- $(a_in-c_im)/2^{s_i}=x_i$
- $(b_in-d_im)/2^{s_i}=-y_i$

## 目標

$\gcd(n,m)=1$ である事を前提に、 $i$ 回の操作後に $x_i=1$ もしくは $y_i=1$ を満たすように式変形が出来たとすると、

$m$ が奇数の場合、モジュラ逆数 $n^{-1}\mod m$ を以下のように導出します。

- $\displaystyle x_i=1\quad\Longrightarrow\quad a_in\equiv 2^{s_i}\quad(\operatorname{mod}m) \quad\Longrightarrow\quad n^{-1}\equiv a_i/2^{s_i}\quad(\operatorname{mod}m)$
- $\displaystyle y_i=1\quad\Longrightarrow\quad b_in\equiv -2^{s_i}\quad(\operatorname{mod}m) \quad\Longrightarrow\quad n^{-1}\equiv -b_i/2^{s_i}\quad(\operatorname{mod}m)$

$m$ が偶数の場合、 $2^{s_i}$ の値が常に $1$ ( つまり、 $s_i=0$ ) となるように変形を行い、 モジュラ逆数 $n^{-1}\mod m$ を以下のように導出します。

- $\displaystyle x_i=1\quad\Longrightarrow\quad a_in\equiv 1\quad(\operatorname{mod}m) \quad\Longrightarrow\quad n^{-1}\equiv a_i\quad(\operatorname{mod}m)$
- $\displaystyle y_i=1\quad\Longrightarrow\quad b_in\equiv -1\quad(\operatorname{mod}m) \quad\Longrightarrow\quad n^{-1}\equiv -b_i\quad(\operatorname{mod}m)$

## 制約

- $\{n,m,a_i,b_i,c_i,d_i,s_i,t_i,x_i,y_i\} は非負整数$
- $\gcd(n,m)=1$

## 初期値

$\displaystyle\begin{pmatrix}s_0&a_0&c_0&x_0\\&b_0&d_0&y_0\end{pmatrix}=\begin{pmatrix}0&1&0&n\\&0&1&m\end{pmatrix}$

- $(a_0n-c_0m)/2^{s_0}=x_0=(1\cdot n-0\cdot m)/1=n$
- $(b_0n-d_0m)/2^{s_0}=-y_0=(0\cdot n-1\cdot m)/1=-m$

## 式変形パターン

A. $x_i \ge t_iy_i$ : 1行目の式に2行目の式の $t_i$ 倍を加える操作

$\gcd(n,m)=\gcd(x_i,y_i)=\gcd(x_i-t_iy_i,y_i)=\gcd(x_{i+1},y_{i+1})$

- $(a_{i+1}n-c_{i+1}m)/2^{s_{i+1}}=x_{i+1}=((a_i+b_it_i)n-(c_i+d_it_i)m)/2^{s_i}=x_i-t_iy_i$
- $(b_{i+1}n-d_{i+1}m)/2^{s_{i+1}}=-y_{i+1}=(b_in-d_im)/2^{s_i}=-y_i$

$$
\begin{pmatrix}s_{i+1}&a_{i+1}&c_{i+1}&x_{i+1}\\&b_{i+1}&d_{i+1}&y_{i+1}\end{pmatrix}=
\begin{pmatrix}s_i&a_i+b_it_i&c_i+d_it_i&x_i-t_iy_i\\&b_i&d_i&y_i\end{pmatrix}
$$

B. $t_ix_i \le y_i$ : 2行目の式に1行目の式の $t_i$ 倍を加える操作

$\gcd(n,m)=\gcd(x_i,y_i)=\gcd(x_i,y_i-t_ix_i)=\gcd(x_{i+1},y_{i+1})$

- $(a_{i+1}n-c_{i+1}m)/2^{s_{i+1}}=x_{i+1}=(a_in-c_im)/2^{s_i}=x_i$
- $(b_{i+1}n-d_{i+1}m)/2^{s_{i+1}}=-y_{i+1}=((a_it_i+b_i)n-(c_it_i+d_i)m)/2^{s_i}=t_ix_i-y_i$

$$
\begin{pmatrix}s_{i+1}&a_{i+1}&c_{i+1}&x_{i+1}\\&b_{i+1}&d_{i+1}&y_{i+1}\end{pmatrix}=
\begin{pmatrix}s_i&a_i&c_i&x_i\\&a_it_i+b_i&c_it_i+d_i&y_i-t_ix_i\end{pmatrix}
$$

C. $m\equiv 1\ (\operatorname{mod}2)$ かつ $x_i$ が $2^{t_i}$ の倍数 : $x_i$の値を $1/2^{t_i}$ する操作

$x_i$ が $2^{t_i}$ の倍数、 $y_i$ が奇数であれば $\gcd(n,m)=\gcd(x_i,y_i)=\gcd(x_i/2^{t_i},y_i)=\gcd(x_{i+1},y_{i+1})$

- $(a_{i+1}n-c_{i+1}m)/2^{s_{i+1}}=x_{i+1}=(a_in-c_im)/2^{s_i+t_i}=x_i/2^{t_i}$
- $(b_{i+1}n-d_{i+1}m)/2^{s_{i+1}}=-y_{i+1}=(2^{t_i}b_in-2^{t_i}d_im)/2^{s_i+t_i}=-y_i$

$$
\begin{pmatrix}s_{i+1}&a_{i+1}&c_{i+1}&x_{i+1}\\&b_{i+1}&d_{i+1}&y_{i+1}\end{pmatrix}=
\begin{pmatrix}s_i+t_i&a_i&c_i&x_i/2^{t_i}\\&2^{t_i}b_i&2^{t_i}d_i&y_i\end{pmatrix}
$$

D. $m\equiv 1\ (\operatorname{mod}2)$ かつ $y_i$ が $2^{t_i}$ の倍数 : $y_i$の値を $1/2^{t_i}$ する操作

$x_i$ が奇数、 $y_i$ が $2^{t_i}$ の倍数であれば $\gcd(n,m)=\gcd(x_i,y_i)=\gcd(x_i,y_i/2^{t_i})=\gcd(x_{i+1},y_{i+1})$

- $(a_{i+1}n-c_{i+1}m)/2^{s_{i+1}}=x_{i+1}=(2^{t_i}a_in-2^{t_i}c_im)/2^{s_i+t_i}=x_i$
- $(b_{i+1}n-d_{i+1}m)/2^{s_{i+1}}=-y_{i+1}=(b_in-d_im)/2^{s_i+t_i}=-y_i/2^{t_i}$

$$
\begin{pmatrix}s_{i+1}&a_{i+1}&c_{i+1}&x_{i+1}\\&b_{i+1}&d_{i+1}&y_{i+1}\end{pmatrix}=
\begin{pmatrix}s_i+t_i&2^{t_i}a_i&2^{t_i}c_i&x_i\\&b_i&d_i&y_i/2^{t_i}\end{pmatrix}
$$

## 式変形手順の例 ( 一般の $m$ について、 A.B. のみの式変形を使う場合 )

A.B. のみの式変形を使う事で、 $\lbrace 2^0, 2^{-1}, 2^{-2}, ..., 2^{-63}\rbrace$ の $\operatorname{mod}m$ 事前計算を省いた手順。

1. $(a,b,x,y)=(1,0,n,m)$ と初期化 ( $c,d,s$ は使わない )
2. $x=0$ であれば異常終了 : $\gcd(n,m)=\gcd(x,y)=y\ge 2$ のためモジュラ逆数が存在しない
3. $x=1$ であれば正常終了 : $\displaystyle n^{-1}\equiv a\quad(\operatorname{mod}m)$ を計算して出力
4. $x\le y$ であれば $t=\lfloor y/x\rfloor$ として 式変形B. を実行
5. $y=0$ であれば異常終了 : $\gcd(n,m)=\gcd(x,y)=x\ge 2$ のためモジュラ逆数が存在しない
6. $y=1$ であれば正常終了 : $\displaystyle n^{-1}\equiv -b\quad(\operatorname{mod}m)$ を計算して出力
7. $x\ge y$ であれば $t=\lfloor x/y\rfloor$ として 式変形A. を実行
8. 手順 2. に戻る

サンプルコード (Rust):

```rust
// 拡張ユークリッド互除法によるモジュラ逆数 $n^{-1} \mod m$
// 符号なし整数型を用いたバリエーション
pub fn modinv_plain(n: u32, m: u32) -> Option<u32> {
    let (mut a, mut b, mut x, mut y) = (1, 0, n, m);
    if m == 1 {
        return Some(0);
    }
    loop {
        match x {
            // $\gcd(n,m) = \gcd(x,y) = y \ge 2$
            0 => return None,
            // $n^{-1} \equiv a$
            1 => return Some(a),
            _ => {}
        }
        b += a * (y / x);
        y %= x;
        match y {
            // $\gcd(n,m) = \gcd(x,y) = x \ge 2$
            0 => return None,
            // $n^{-1} \equiv -b$
            1 => return Some(m - b),
            _ => {}
        }
        a += b * (x / y);
        x %= y;
    }
}
```

## 式変形手順の具体例 ( $m$ が $3$ 以上の奇数 の場合、バイナリGCD法の複合 )

$x,y$ の比が大きい時のみ剰余算を用いた A.B. の式変形を使い、比が小さい時は 加減算を用いた A.B. の式変形 と C.D. の式変形 (バイナリGCD法) を使う事で、それぞれの手順の性能の美味しいところ取りを狙ったもの。

https://en.wikipedia.org/wiki/Binary_GCD_algorithm

1. 例えば $0\le n\lt m\lt 2^{32}$ の場合、 $\lbrace 2^0, 2^{-1}, 2^{-2}, ..., 2^{-63}\rbrace$ の $\operatorname{mod}m$ をそれぞれ事前計算しておく ( 計算方法は後述の補題参照、 $\lfloor\log_2 mn\rfloor$ が取りうる値くらいの数だけ用意 )
2. $(a,b,s,x,y)=(1,0,0,n,m)$ と初期化 ( $c,d$ は使わない )
3. $x=0$ であれば異常終了 : $\gcd(n,m)=\gcd(x,y)=y\ge 3$ のためモジュラ逆数が存在しない
4. $x\ge 2$ の偶数であれば $x$ を奇数にするような 自然数$t$ にて 式変形C. を実行
5. $x=1$ であれば正常終了 : $\displaystyle n^{-1}\equiv 2^{-s}a\quad(\operatorname{mod}m)$ を計算して出力
6. $x=y$ であれば異常終了 : $\gcd(n,m)=\gcd(x,y)=x=y\ge 3$ のためモジュラ逆数が存在しない
7. $x\lt 8y$ であれば $t=\lfloor y/x\rfloor$ として 式変形B. を実行
8. $x\lt y$ であれば $t=1$ として式変形B. を実行
9. $y=0$ であれば異常終了 : $\gcd(n,m)=\gcd(x,y)=x\ge 3$ のためモジュラ逆数が存在しない
10. $y\ge 2$ の偶数であれば $y$ を奇数にするような 自然数$t$ にて 式変形D. を実行
11. $y=1$ であれば正常終了 : $\displaystyle n^{-1}\equiv -2^{-t}b\quad(\operatorname{mod}m)$ を計算して出力
12. $8x\gt y$ であれば $t=\lfloor x/y\rfloor$ として 式変形A. を実行
13. $x\gt y$ であれば $t=1$ として式変形A. を実行
14. 手順 3. に戻る

サンプルコード (Rust): 記事末尾のサンプルコード `fn inv_mx` 関数の部分を参照。

## 補題: $m$ が $3$ 以上の奇数の場合の $2$ の累乗数に対するモジュラ逆数の事前計算

$m$ が奇数 : $m\equiv 1\ (\operatorname{mod}2)$ の時、

- $1/2^0\ \operatorname{mod}m=1$
- $\displaystyle 1/2^{q+1}\ \operatorname{mod}m=\begin{cases}(1/2^{q}\ \operatorname{mod}m)/2&((1/2^{q}\ \operatorname{mod}m)\equiv 0\ \operatorname{mod}2)\\ ((1/2^{q}\ \operatorname{mod}m)+m)/2&((1/2^{q}\ \operatorname{mod}m)\equiv 1\ \operatorname{mod}2)\end{cases}$

例:

$$
\begin{align*}
1 \times 2^{0} &\equiv 1 \mod 998244353\\
(1+998244353)/2 \times 2^{1} \equiv 499122177 \times 2^{1} &\equiv 1 \mod 998244353\\
(499122177+998244353)/2 \times 2^{2} \equiv 748683265 \times 2^{2} &\equiv 1 \mod 998244353\\
(748683265+998244353)/2 \times 2^{3} \equiv 873463809 \times 2^{3} &\equiv 1 \mod 998244353\\
(873463809+998244353)/2 \times 2^{4} \equiv 935854081 \times 2^{4} &\equiv 1 \mod 998244353\\
(935854081+998244353)/2 \times 2^{5} \equiv 967049217 \times 2^{5} &\equiv 1 \mod 998244353\\
(967049217+998244353)/2 \times 2^{6} \equiv 982646785 \times 2^{6} &\equiv 1 \mod 998244353\\
&\vdots\\
(436731897+998244353)/2 \times 2^{28} \equiv 717488125 \times 2^{28} &\equiv 1 \mod 998244353\\
(717488125+998244353)/2 \times 2^{29} \equiv 857866239 \times 2^{29} &\equiv 1 \mod 998244353\\
(857866239+998244353)/2 \times 2^{30} \equiv 928055296 \times 2^{30} &\equiv 1 \mod 998244353\\
(928055296)/2 \times 2^{31} \equiv 464027648 \times 2^{31} &\equiv 1 \mod 998244353\\
(464027648)/2 \times 2^{32} \equiv 232013824 \times 2^{32} &\equiv 1 \mod 998244353\\
(232013824)/2 \times 2^{33} \equiv 116006912 \times 2^{33} &\equiv 1 \mod 998244353\\
&\vdots\\
\end{align*}
$$

## $m=1$ の場合

1. $nn^{-1} \equiv 0 \equiv 1\ (\operatorname{mod}1)$ より、 $n^{-1}\equiv 0\ (\operatorname{mod}1)$ を出力

## $m=0$ の場合

1. 未定義動作とする。異常終了しても良いし、 $m=2^{32}$ や $m=2^{64}$ などとみなして処理しても良いし、鼻から悪魔が出ても良い。

## ニュートン・ラフソン法、中国剰余定理 を組み合わせる (任意の $m$ の場合、主に $m$ が $2$ 以上の偶数 の場合)

$m$ が $2$ 以上の偶数の場合に、 モジュラ逆数 $n^{-1}\ \operatorname{mod}m$ を求める手順の例。ここでは前述の 「ニュートン・ラフソン法」、「変則的拡張ユークリッド互除法・式変形手順の具体例 ( $m$ が $3$ 以上の奇数 の場合 )」 に後述の 「中国剰余定理 (Garner のアルゴリズム)」 を組み合わせ、さらに変則的な拡張ユークリッド互除法を構成してみます。

1. $m=m_xm_y,\ m_x=d,\ m_y=2^r$ とする。 ここで、 $r,d$ は $m=2^r\times d$ の形に分解した時の値 ( $r$ は非負整数、 $d$ は正の奇数 ) 。 $m_x$ は正の奇数、 $m_y$ は $2$ の累乗数 ( ただし $m$ が奇数であれば、 $m_y=1$ ) であり、 $\gcd(m_x, m_y)=1$ を満たす。また、 $\gcd(n,m)=1$ であれば、 $\gcd(n,m_x)=1,\ \gcd(n,m_y)=1$ を満たす。
2. $m_x\overline{m_x}\equiv 1\quad(\operatorname{mod}m_y)$ を満たす、 $\operatorname{mod}m_y$ における $m_x$ のモジュラ逆数 $\overline{m_x}$ の値を、前述の「ニュートン・ラフソン法」で求める。 $m_y=1$ の時は、 $\overline{m_x}=0$ とする。
3. $\operatorname{mod}m_x$ における $n$ のモジュラ逆数 $x=n^{-1}\ \operatorname{mod}m_x$ の値を、前述の「変則的拡張ユークリッドの互除法・式変形手順の具体例 ( $m$ が $3$ 以上の奇数 の場合 )」で求める。 $m_x=1$ の時は、 $x=0$ とする。 $\gcd(n,m_x)\ne 1$ であれば、 $\gcd(n,m)\ne 1$ であるので、モジュラ逆数 $n^{-1}\ \operatorname{mod}m$ の値が存在しないとして異常終了する。
4. $\operatorname{mod}m_y$ における $n$ のモジュラ逆数 $y=n^{-1}\ \operatorname{mod}m_y$ の値を、前述の「ニュートン・ラフソン法」で求める。 $m_y=1$ の時は、 $y=0$ とする。 $\gcd(n,m_y)\ne 1$ であれば、 $\gcd(n,m)\ne 1$ であるので、モジュラ逆数 $n^{-1}\ \operatorname{mod}m$ の値が存在しないとして異常終了する。
5. $m=m_xm_y,\ \gcd(m_x,m_y)=1$ であるので、 $z\equiv x\ (\operatorname{mod}m_x),\ z\equiv y\ (\operatorname{mod}m_y)$ をそれぞれ満たすような 整数 $0\le z\lt m$ の値を、後述の「中国剰余定理 (Garner のアルゴリズム)」を用いて、 $z=(x-m_x((\overline{m_x}(x-y))\operatorname{mod}m_y))\operatorname{mod}m$ もしくは $z=(x+m_x((\overline{m_x}(y-x))\operatorname{mod}m_y))\operatorname{mod}m$ のように求める。
6. $m=m_xm_y,\ \gcd(m_x,m_y)=1, nz\equiv 1\ (\operatorname{mod}m_x),\ nz\equiv 1\ (\operatorname{mod}m_y)$ より $nz\equiv 1\ (\operatorname{mod}m)$ であるので、 $z$ を モジュラ逆数 $n^{-1}\ \operatorname{mod}m$ の値として正常終了する。

## 補題: 中国剰余定理 (Garner のアルゴリズム)

https://ja.wikipedia.org/wiki/%E4%B8%AD%E5%9B%BD%E3%81%AE%E5%89%B0%E4%BD%99%E5%AE%9A%E7%90%86

https://qiita.com/drken/items/ae02240cd1f8edfc86fd#2-2-garner-%E3%81%AE%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0

$\gcd(M_x,M_y)=1$ の時、 $z\equiv x\ (\operatorname{mod}M_x),\ z\equiv y\ (\operatorname{mod}M_y)$ が成り立つような整数 $0\le z\lt M_xM_y$ の値を求めるアルゴリズム。

- $x,y\in\text{整数}\mathbb{Z}$
- $M_x,M_y,\overline{M_x}\in\text{自然数}\mathbb{N}$
- $\gcd(M_x,M_y)=1$
- $\overline{M_x}M_x\equiv 1\quad(\operatorname{mod}M_y)$
- $t=M_x((\overline{M_x}(x-y))\operatorname{mod}M_y)$
- $z=\begin{cases}x-t&(x \ge t)\\x-t+M_xM_y&(x \lt t)\\\end{cases}$

とすると

- $0\le t\le M_x(M_y-1)$
- $x-M_x(M_y-1)\le x-t\le x+M_xM_y$
- $z\equiv x\quad(\operatorname{mod}M_x)\quad\because\quad t\equiv 0\quad(\operatorname{mod}M_x)$
- $z\equiv y\quad(\operatorname{mod}M_y)\quad\because\quad t\equiv x-y\quad(\operatorname{mod}M_y)$
- 特に、 $0\le x\lt M_xM_y$ であれば $0\le z\lt M_xM_y$

## サンプルコード (Rust)

拡張ユークリッド互除法 + バイナリGCD法 + ニュートン・ラフソン法 + Garnerのアルゴリズム を組み合わせてモジュラ逆元を計算するコードの例です。

```rust
// Modulus値分類
pub enum ModIntDynType {
    Zero,
    Odd,
    P0998244353,
    P1000000007,
    Power2,
    Even,
}
// モジュラ演算用定数構造体
pub struct ModIntDyn {
    ty: ModIntDynType,
    m: u32,
    tz: u32,
    mx: u32,
    mimx: u32,
    mym: u32,
    dim: u64,
    dimx: u64,
    ib: [u32; 64],
}
impl ModIntDyn {
    pub fn new_modulus(m: u32) -> Self {
        if m == 0 {
            // $m = 0$ であれば $m = 2^{32}$ とみなして処理
            return Self {
                m,
                ty: ModIntDynType::Zero,
                tz: 32,
                mx: 1,
                mimx: 0,
                mym: !0u32,
                dim: 1,
                dimx: 0,
                ib: [1u32; 64],
            };
        }
        // tz = tzcnt(m), my = 2^tz
        let tz = m.trailing_zeros();
        // mx = modulus / my
        let mx = m >> tz;
        // mx * imx % 2^{32} = 1 => mx * imx % my = 1 % my
        let mimx = {
            let mut t = mx.wrapping_mul(3) ^ 2;
            let mut u = 1u32.wrapping_sub(t.wrapping_mul(mx));
            t = t.wrapping_mul(u.wrapping_add(1));
            u = u.wrapping_mul(u);
            t = t.wrapping_mul(u.wrapping_add(1));
            u = u.wrapping_mul(u);
            t = t.wrapping_mul(u.wrapping_add(1));
            t
        };
        // mym = my - 1, x % my == x & mym
        let mym = u32::from(tz != 0).wrapping_neg() >> ((32 - tz) & 31);
        // `dim` $\ = \lceil (2^{64} / m) \rceil = \lfloor (2^{64}-1) / m \rfloor + 1$
        let dim = ((!0u64) / u64::from(m)).wrapping_add(1);
        // `dimx` $\ = \lceil (2^{64} / m_x) \rceil = \lfloor (2^{64}-1) / m_x \rfloor + 1$
        let dimx = ((!0u64) / u64::from(mx)).wrapping_add(1);
        // ib: $\{2^0,2^{-1},\dots,2^{-63}\} \mod m_x$ を事前計算
        let mut ib = [0u32; 64];
        let mut t = 1;
        for e in ib.iter_mut() {
            *e = t;
            t = t / 2 + ((t % 2).wrapping_neg() & (mx / 2 + 1));
        }
        // オブジェクト生成
        Self {
            m,
            ty: {
                if m == 998_244_353 {
                    ModIntDynType::P0998244353
                } else if m == 1_000_000_007 {
                    ModIntDynType::P1000000007
                } else if tz == 0 {
                    ModIntDynType::Odd
                } else if mx == 1 {
                    ModIntDynType::Power2
                } else {
                    ModIntDynType::Even
                }
            },
            tz,
            mx,
            mimx,
            mym,
            dim,
            dimx,
            ib,
        }
    }
    // $a \mod m$ (normalize)
    #[inline]
    pub fn norm(&self, a: u32) -> u32 {
        match self.ty {
            // $m = 2^{32}$
            ModIntDynType::Zero => a,
            // $m = m_y$
            ModIntDynType::Power2 => a & self.mym,
            // others
            _ => a % self.m,
        }
    }
    // $a + b \mod m$
    #[inline]
    pub fn add(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.m);
        debug_assert!(b < self.m);
        let v = u64::from(a) + u64::from(b);
        match self.ty {
            // $m = 2^{32}
            ModIntDynType::Zero => v as u32,
            // others
            _ => match v.overflowing_sub(u64::from(self.m)) {
                (_, true) => v as u32,
                (w, false) => w as u32,
            },
        }
    }
    // $a - b \mod m$
    #[inline]
    pub fn sub(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.m);
        debug_assert!(b < self.m);
        let (v, f) = a.overflowing_sub(b);
        match self.ty {
            // $m = 2^{32}
            ModIntDynType::Zero => v,
            // others
            _ => v.wrapping_add(u32::from(f).wrapping_neg() & self.m),
        }
    }
    // $-a \mod m$
    #[inline]
    pub fn neg(&self, a: u32) -> u32 {
        debug_assert!(a < self.m);
        self.sub(0, a)
    }
    // $a * b \mod m$
    #[inline]
    pub fn mul(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.m);
        debug_assert!(b < self.m);
        let z = u64::from(a) * u64::from(b);
        match self.ty {
            // $m = m_y$
            ModIntDynType::Power2 => (z as u32) & self.mym,
            // $m = 998244353$
            ModIntDynType::P0998244353 => {
                (z - (((u128::from(z) * 0x89AE_4087_5DE0_CC3F) >> 93) as u64) * 998_244_353) as u32
            }
            // $m = 1000000007$
            ModIntDynType::P1000000007 => {
                (z - (((u128::from(z) * 0x8970_5F31_12A2_8FE5) >> 93) as u64) * 1_000_000_007)
                    as u32
            }
            // $m = 2^{32}$
            ModIntDynType::Zero => z as u32,
            // others
            ModIntDynType::Odd | ModIntDynType::Even => {
                // Barrett Reduction
                let x = ((u128::from(z) * u128::from(self.dim)) >> 64) as u64;
                let (v, f) = z.overflowing_sub(x.wrapping_mul(u64::from(self.m)));
                (v as u32).wrapping_add(u32::from(f).wrapping_neg() & self.m)
            }
        }
    }
    // $a * b \mod m_x$
    #[inline]
    fn mulmx(&self, a: u32, b: u32) -> u32 {
        debug_assert!(a < self.mx);
        debug_assert!(b < self.mx);
        let z = u64::from(a) * u64::from(b);
        match self.ty {
            // $m = 998244353$
            ModIntDynType::P0998244353 => {
                (z - (((u128::from(z) * 0x89AE_4087_5DE0_CC3F) >> 93) as u64) * 998_244_353) as u32
            }
            // $m = 1000000007$
            ModIntDynType::P1000000007 => {
                (z - (((u128::from(z) * 0x8970_5F31_12A2_8FE5) >> 93) as u64) * 1_000_000_007)
                    as u32
            }
            // $m_x = 0$
            ModIntDynType::Power2 | ModIntDynType::Zero => 0,
            // others
            ModIntDynType::Odd | ModIntDynType::Even => {
                // Barrett Reduction
                let x = ((u128::from(z) * u128::from(self.dimx)) >> 64) as u64;
                let (v, f) = z.overflowing_sub(x.wrapping_mul(u64::from(self.mx)));
                (v as u32).wrapping_add(u32::from(f).wrapping_neg() & self.mx)
            }
        }
    }
    // $a^b \mod m$
    #[inline]
    pub fn pow(&self, a: u32, mut b: u32) -> u32 {
        debug_assert!(a < self.m);
        if b == 0 {
            return (self.m > 1) as u32;
        }
        let mut t = a;
        while b % 2 == 0 {
            t = self.mul(t, t);
            b /= 2;
        }
        let mut r = t;
        b /= 2;
        while b > 0 {
            t = self.mul(t, t);
            if b % 2 != 0 {
                r = self.mul(r, t);
            }
            b /= 2;
        }
        r
    }
    // $a / b \mod m$
    #[inline]
    pub fn div(&self, a: u32, b: u32) -> Option<u32> {
        debug_assert!(a < self.m);
        debug_assert!(b < self.m);
        self.inv(b).map(|ib| self.mul(a, ib))
    }
    // モジュラ逆数の計算(非負整数による拡張ユークリッド互除法の実装例)
    #[inline]
    pub fn inv_plain(&self, n: u32) -> Option<u32> {
        let m = self.m;
        let (mut a, mut b, mut x, mut y) = (1, 0, n, m);
        if m == 1 {
            return Some(0);
        }
        loop {
            match x {
                // $\gcd(n,m) = \gcd(x,y) = y \ge 2$
                0 => return None,
                // $n^{-1} \equiv a$
                1 => return Some(a),
                _ => {}
            }
            b += a * (y / x);
            y %= x;
            match y {
                // $\gcd(n,m) = \gcd(x,y) = x \ge 2$
                0 => return None,
                // $n^{-1} \equiv -b$
                1 => return Some(m - b),
                _ => {}
            }
            a += b * (x / y);
            x %= y;
        }
    }
    // inv の子関数: 拡張ユークリッド互除法 + バイナリGCD法 による モジュラ逆数の計算 ($\mod m_x$) ($m_x$ は 正の奇数)
    #[inline]
    fn inv_mx(&self, n: u32) -> Option<u32> {
        debug_assert_ne!(self.mx % 2, 0);
        if self.mx % 2 == 0 {
            unsafe { core::hint::unreachable_unchecked() }
        }
        if self.mx == 1 {
            return Some(0);
        }
        let (mut a, mut b, mut x, mut y) = (1, 0, n, self.mx);
        if x == 0 {
            // $\gcd(n,m_x) = \gcd(x,y) = m_x \ge 3$
            return None;
        }
        // 式変形C: $t = \operatorname{tzcnt}(x)
        let mut s = x.trailing_zeros();
        x >>= s;
        loop {
            if x == 1 {
                // $n^{-1} \equiv a / 2^s \mod m_x$
                return Some(self.mulmx(a, self.ib[s as usize]));
            }
            if x == y {
                // $\gcd(n,m_x) = \gcd(x,y) = x = y \ge 3$
                return None;
            }
            if x <= (y >> 3) {
                // 式変形B: $t = \lfloor y / x \rfloor$
                if x == 0 {
                    unsafe { core::hint::unreachable_unchecked() }
                }
                b += a * (y / x);
                y %= x;
                if y == 0 {
                    // $\gcd(n,m_x) = \gcd(x,y) = x \ge 3$
                    return None;
                }
                // 式変形D: $t = \operatorname{tzcnt}(y)$
                let t = y.trailing_zeros();
                y >>= t;
                a <<= t;
                s += t;
            } else {
                while x < y {
                    // 式変形B: $t = 1$
                    b += a;
                    y -= x;
                    // 式変形D: $t = \operatorname{tzcnt}(y)$
                    let t = y.trailing_zeros();
                    y >>= t;
                    a <<= t;
                    s += t;
                }
            }
            if y == 1 {
                // $n^{-1} \equiv -b / 2^s \mod m_x$
                return Some(self.mulmx(self.mx - b, self.ib[s as usize]));
            }
            if (x >> 3) >= y {
                // 式変形A: $t = \lfloor x / y \rfloor$
                if y == 0 {
                    unsafe { core::hint::unreachable_unchecked() }
                }
                a += b * (x / y);
                x %= y;
                if x == 0 {
                    // $\gcd(n,m_x) = \gcd(x,y) = y \ge 3$
                    return None;
                }
                // 式変形C: $t = \operatorname{tzcnt}(x)$
                let t = x.trailing_zeros();
                x >>= t;
                b <<= t;
                s += t;
            } else {
                while x > y {
                    // 式変形A: $t = 1$
                    a += b;
                    x -= y;
                    // 式変形C: $t = \operatorname{tzcnt}(x)$
                    let t = x.trailing_zeros();
                    x >>= t;
                    b <<= t;
                    s += t;
                }
            }
        }
    }
    // inv の子関数: ニュートン・ラフソン法 による モジュラ逆数 の計算 ($\mod m_y$) ($m_y$ は 1 もしくは 2 の累乗数)
    #[inline]
    fn inv_my(&self, n: u32) -> Option<u32> {
        if self.mym == 0 {
            // 0bit, $m_y == 1$
            return Some(0);
        }
        if n % 2 == 0 {
            // $\gcd(n,m_y) \ge 2$
            return None;
        }
        Some(
            match self.tz {
                // 21bit ~ 32bit
                21..=0xffff_ffff => {
                    let mut t = n.wrapping_mul(3) ^ 2;
                    let mut u = 1u32.wrapping_sub(t.wrapping_mul(n));
                    t = t.wrapping_mul(u.wrapping_add(1));
                    u = u.wrapping_mul(u);
                    t = t.wrapping_mul(u.wrapping_add(1));
                    u = u.wrapping_mul(u);
                    t = t.wrapping_mul(u.wrapping_add(1));
                    t
                }
                // 11bit ~ 20bit
                11..=20 => {
                    let mut t = n.wrapping_mul(3) ^ 2;
                    let mut u = 1u32.wrapping_sub(t.wrapping_mul(n));
                    t = t.wrapping_mul(u.wrapping_add(1));
                    u = u.wrapping_mul(u);
                    t = t.wrapping_mul(u.wrapping_add(1));
                    t
                }
                // 6bit ~ 10bit
                6..=10 => {
                    let t = n.wrapping_mul(3) ^ 2;
                    t.wrapping_mul(2u32.wrapping_sub(n.wrapping_mul(t)))
                }
                // 4bit ~ 5bit
                4..=5 => n.wrapping_mul(3) ^ 2,
                // 1bit ~ 3bit
                0..=3 => n,
            } & self.mym,
        )
    }
    // モジュラ逆数の計算(拡張ユークリッドの互除法、バイナリGCD法、ニュートン・ラフソン法、Garnerのアルゴリズムのハイブリッド)
    #[inline]
    pub fn inv(&self, n: u32) -> Option<u32> {
        match self.ty {
            ModIntDynType::Power2 => self.inv_my(n),
            ModIntDynType::P0998244353 | ModIntDynType::P1000000007 | ModIntDynType::Odd => {
                self.inv_mx(n)
            }
            ModIntDynType::Even => {
                if let Some(y) = self.inv_my(n) {
                    // 拡張ユークリッドの互除法+バイナリGCD法 にて $n^{-1} \mod m_x$ を計算
                    if let Some(x) = self.inv_mx(n) {
                        // 中国剰余定理・Garner のアルゴリズム による結果の合成
                        let (t, f) = x.overflowing_sub(
                            self.mx * (self.mimx.wrapping_mul(x.wrapping_sub(y)) & self.mym),
                        );
                        return Some(t.wrapping_add(u32::from(f).wrapping_neg() & self.m));
                    }
                }
                None
            }
            ModIntDynType::Zero => {
                if n % 2 == 0 {
                    return None;
                }
                let mut t = n.wrapping_mul(3) ^ 2;
                let mut u = 1u32.wrapping_sub(t.wrapping_mul(n));
                t = t.wrapping_mul(u.wrapping_add(1));
                u = u.wrapping_mul(u);
                t = t.wrapping_mul(u.wrapping_add(1));
                u = u.wrapping_mul(u);
                t = t.wrapping_mul(u.wrapping_add(1));
                Some(t)
            }
        }
    }
}

// 1 / values[0] / values[1] / ... (mod m)
pub fn moddyn_div_test(modulus: u32, values: &[u32]) -> u32 {
    let moddyn = ModIntDyn::new_modulus(modulus);
    values.iter().fold(1u32, |p, &c| moddyn.div(p, c).unwrap())
}
```
