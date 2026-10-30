# Rzut ukośny --- model bez oporu i z oporem powietrza

*Raport obliczeniowy: wyprowadzenie równań, symulacja numeryczna i
interpretacja wykresów.*

**Temat:** ruch ciała wyrzuconego pod kątem do poziomu.  
**Parametry początkowe:** $v_0 = 36{,}0\,\mathrm{m/s}$, $\alpha = 79{,}0^\circ$.  
**Metoda:** rozwiązanie analityczne dla próżni oraz numeryczne całkowanie równań ruchu z oporem kwadratowym.  
**Charakter danych:** wyniki modelowe, nie pomiary doświadczalne.

Rzut ukośny jest przykładem ruchu, który można rozłożyć na dwa
niezależne kierunki tylko wtedy, gdy pomijamy opór powietrza. W modelu
idealnym składowa pozioma prędkości pozostaje stała, a pionowa zmienia
się pod wpływem grawitacji. Po uwzględnieniu oporu obie składowe są
sprzężone przez wartość całkowitej prędkości.

## Spis treści

1.  [Cel i założenia modelu](#1-cel-i-założenia-modelu)
2.  [Rzut ukośny bez oporu powietrza](#2-rzut-ukośny-bez-oporu-powietrza)
3.  [Przykład liczbowy --- model idealny](#3-przykład-liczbowy--model-idealny)
4.  [Rzut ukośny z oporem powietrza](#4-rzut-ukośny-z-oporem-powietrza)
5.  [Porównanie wyników i interpretacja wykresów](#5-porównanie-wyników-i-interpretacja-wykresów)
6.  [Ograniczenia modelu](#6-ograniczenia-modelu)
7.  [Kod i odtwarzanie wykresów](#7-kod-i-odtwarzanie-wykresów)
8.  [Wnioski](#8-wnioski)
9.  [Źródła i dalsza lektura](#9-źródła-i-dalsza-lektura)

------------------------------------------------------------------------

## 1. Cel i założenia modelu

Celem raportu jest opisanie toru, czasu lotu, maksymalnej wysokości,
zasięgu, prędkości i energii ciała wyrzuconego z poziomu gruntu pod
kątem $\alpha$. Rozpatrujemy dwa przypadki:

-   **model idealny** --- brak oporu powietrza, stałe przyspieszenie
    grawitacyjne;
-   **model z oporem** --- siła oporu aerodynamicznego zależy od
    kwadratu szybkości i jest skierowana przeciwnie do wektora
    prędkości.

### 1.1. Przyjęte dane

| Wielkość | Symbol | Wartość |
| --- | --- | --- |
| Prędkość początkowa | $v_0$ | $36{,}0\,\mathrm{m/s}$ |
| Kąt wyrzutu | $\alpha$ | $79{,}0^\circ$ |
| Przyspieszenie ziemskie | $g$ | $9{,}81\,\mathrm{m/s^2}$ |
| Masa ciała w modelu oporu | $m$ | $0{,}145\,\mathrm{kg}$ |
| Gęstość powietrza | $\rho$ | $1{,}225\,\mathrm{kg/m^3}$ |
| Współczynnik oporu | $C_D$ | $0{,}47$ |
| Średnica ciała | $d$ | $0{,}073\,\mathrm{m}$ |

Wartości $m,\rho,C_D,d$ są dodatkowymi założeniami modelu
potrzebnymi do obliczenia oporu. Nie wynikają z podanych parametrów
indeksu. Przyjęto przykładowe parametry niewielkiej kulistej piłki;
rzeczywisty współczynnik oporu zależy m.in. od kształtu, powierzchni i
liczby Reynoldsa.

### 1.2. Układ współrzędnych i warunki początkowe

Początek układu współrzędnych znajduje się w punkcie wyrzutu. Oś $x$
jest pozioma i skierowana w prawo, a oś $y$ pionowa i skierowana w górę.
Poziom gruntu ma równanie $y=0$.

W chwili $t=0$:

$$
 x(0)=0,\qquad y(0)=0.
$$

Składowe prędkości początkowej wynoszą:

$$
 v_{x0}=v_0\cos\alpha,\qquad v_{y0}=v_0\sin\alpha.
$$

Po podstawieniu danych:

$$
 v_{x0}=36\cos79^\circ\approx 6{,}87\,\mathrm{m/s},
$$

$$
 v_{y0}=36\sin79^\circ\approx 35{,}34\,\mathrm{m/s}.
$$

### 1.3. Symbole i jednostki SI

| Symbol | Znaczenie | Jednostka |
| --- | --- | --- |
| $(x,y)$ | współrzędne położenia | $\mathrm{m}$ |
| $t$ | czas | $\mathrm{s}$ |
| $(v_x,v_y,v)$ | składowe i wartość prędkości | $\mathrm{m/s}$ |
| $(a_x,a_y)$ | składowe przyspieszenia | $\mathrm{m/s^2}$ |
| $m$ | masa | $\mathrm{kg}$ |
| $g$ | przyspieszenie ziemskie | $\mathrm{m/s^2}$ |
| $\rho$ | gęstość powietrza | $\mathrm{kg/m^3}$ |
| $C_D$ | współczynnik oporu | $1$ |
| $A$ | pole przekroju poprzecznego | $\mathrm{m^2}$ |
| $F_D$ | wartość siły oporu | $\mathrm{N}$ |
| $E$ | energia mechaniczna | $\mathrm{J}$ |

------------------------------------------------------------------------

## 2. Rzut ukośny bez oporu powietrza

### 2.1. Równania ruchu

W modelu idealnym na ciało działa wyłącznie siła ciężkości:

$$
\vec{F}_g = m\vec{g}.
$$

Z drugiej zasady dynamiki Newtona:

$$
 m\vec{a}=m\vec{g}.
$$

Dla osi poziomej i pionowej otrzymujemy:

$$
 m a_x=0,\qquad m a_y=-mg.
$$

Po podzieleniu przez masę:

$$
\boxed{a_x=0,\qquad a_y=-g.}
$$

Ponieważ przyspieszenie jest pochodną prędkości po czasie:

$$
\frac{dv_x}{dt}=0,\qquad \frac{dv_y}{dt}=-g.
$$

Całkując i uwzględniając prędkość początkową:

$$
\boxed{v_x(t)=v_0\cos\alpha,}
$$

$$
\boxed{v_y(t)=v_0\sin\alpha-gt.}
$$

### 2.2. Równania położenia

Prędkość jest pochodną położenia:

$$
\frac{dx}{dt}=v_0\cos\alpha,
$$

$$
\frac{dy}{dt}=v_0\sin\alpha-gt.
$$

Całkujemy pierwsze równanie:

$$
 x(t)=v_0\cos\alpha\cdot t+C_x.
$$

Z warunku $x(0)=0$ wynika $C_x=0$, zatem:

$$
\boxed{x(t)=v_0\cos\alpha\cdot t.}
$$

Analogicznie dla osi pionowej:

$$
 y(t)=v_0\sin\alpha\cdot t-\frac{1}{2}gt^2 + C_y.
$$

Z warunku $y(0)=0$ otrzymujemy $C_y=0$:

$$
\boxed{y(t)=v_0\sin\alpha\cdot t-\frac{1}{2}gt^2.}
$$

### 2.3. Równanie toru

Aby otrzymać równanie toru, eliminujemy czas. Z równania poziomego:

$$
 t=\frac{x}{v_0\cos\alpha}.
$$

Podstawiamy do równania pionowego:

$$
 y=x\tan\alpha - \frac{g x^2}{2v_0^2\cos^2\alpha}.
$$

Jest to równanie paraboli:

$$
\boxed{y(x)=x\tan\alpha-\frac{g}{2v_0^2\cos^2\alpha}x^2.}
$$

W modelu bez oporu tor ruchu jest symetryczną parabolą, jeżeli punkt
lądowania znajduje się na tej samej wysokości co punkt wyrzutu.

### 2.4. Czas wznoszenia i maksymalna wysokość

W najwyższym punkcie toru pionowa składowa prędkości jest równa zeru:

$$
 v_y(t_{\max})=0.
$$

Stąd:

$$
 v_0\sin\alpha - gt_{\max}=0,
$$

$$
\boxed{t_{\max}=\frac{v_0\sin\alpha}{g}.}
$$

Maksymalną wysokość otrzymujemy przez podstawienie tego czasu do $y(t)$:

$$
 H_{\max}=v_0\sin\alpha\frac{v_0\sin\alpha}{g} - \frac{g}{2}\left(\frac{v_0\sin\alpha}{g}\right)^2.
$$

Po uproszczeniu:

$$
\boxed{H_{\max}=\frac{v_0^2\sin^2\alpha}{2g}.}
$$

### 2.5. Czas lotu i zasięg

Ciało wraca na poziom gruntu, gdy $y(t)=0$:

$$
 v_0\sin\alpha\cdot t - \frac{1}{2}gt^2 =0.
$$

Wyłączamy $t$:

$$
 t\left(v_0\sin\alpha - \frac{1}{2}gt\right)=0.
$$

Rozwiązanie $t=0$ oznacza chwilę wyrzutu. Drugie rozwiązanie to czas
lotu:

$$
\boxed{T=\frac{2v_0\sin\alpha}{g}.}
$$

Zasięg poziomy jest równy $x(T)$:

$$
 R=v_0\cos\alpha\cdot T,
$$

czyli:

$$
\boxed{R=\frac{v_0^2\sin(2\alpha)}{g}.}
$$

Wzór ten obowiązuje dla startu i lądowania na tej samej wysokości oraz
przy braku oporu powietrza.

------------------------------------------------------------------------

## 3. Przykład liczbowy --- model idealny

Dane:

$$
 v_0=36\,\mathrm{m/s},\qquad \alpha=79^\circ,\qquad g=9{,}81\,\mathrm{m/s^2}. 
$$

### 3.1. Wyniki obliczeń

| Wielkość | Wzór | Wynik |
| --- | --- | --- |
| Pozioma składowa prędkości początkowej | $v_0\cos\alpha$ | ok. $6{,}87\,\mathrm{m/s}$ |
| Pionowa składowa prędkości początkowej | $v_0\sin\alpha$ | ok. $35{,}34\,\mathrm{m/s}$ |
| Czas wznoszenia | $v_0\sin\alpha/g$ | ok. $3{,}60\,\mathrm{s}$ |
| Czas całego lotu | $2v_0\sin\alpha/g$ | ok. $7{,}20\,\mathrm{s}$ |
| Maksymalna wysokość | $v_0^2\sin^2\alpha/(2g)$ | ok. $63{,}6\,\mathrm{m}$ |
| Zasięg | $v_0^2\sin(2\alpha)/g$ | ok. $49{,}5\,\mathrm{m}$ |

### 3.2. Przebiegi czasowe

Dla przyjętych danych równania ruchu mają postać:

$$
 x(t)=36\cos79^\circ\cdot t,
$$

$$
 y(t)=36\sin79^\circ\cdot t-4{,}905t^2.
$$

Składowe prędkości:

$$
 v_x(t)=36\cos79^\circ,
$$

$$
 v_y(t)=36\sin79^\circ-9{,}81t.
$$

Składowe przyspieszenia są stałe:

$$
 a_x(t)=0,\qquad a_y(t)=-9{,}81\,\mathrm{m/s^2}.
$$

Wykres toru pokazuje parabolę. Wykresy położenia i prędkości pozwalają
odczytać moment osiągnięcia najwyższego punktu. W tej chwili $v_y=0$,
ale $v_x$ nadal ma niezerową wartość.

![Tor ruchu w modelu bez oporu i z oporem powietrza](files/rzut_ukosny/01_tor_ruchu.png)

*Rysunek 1. Porównanie torów ruchu w obu modelach.*

![Położenie poziome i pionowe w funkcji czasu](files/rzut_ukosny/02_polozenie_czas.png)

*Rysunek 2. Zmiana współrzędnych położenia w czasie.*

![Składowe prędkości w funkcji czasu](files/rzut_ukosny/03_predkosci.png)

*Rysunek 3. Składowe prędkości poziomej i pionowej.*

------------------------------------------------------------------------

## 4. Rzut ukośny z oporem powietrza

### 4.1. Model siły oporu

W rzeczywistym ruchu powietrze oddziałuje na poruszające się ciało. Dla
odpowiednio dużych szybkości często stosuje się model oporu
kwadratowego:

$$
\boxed{\vec F_D=-\frac12\rho C_D A|\vec v|\vec v.}
$$

Znak minus oznacza, że siła oporu ma zwrot przeciwny do prędkości.
Wartość siły wynosi:

$$
 F_D=\frac12\rho C_D A v^2.
$$

Pole przekroju kulistego ciała:

$$
 A=\pi r^2 = \pi\left(\frac{d}{2}\right)^2.
$$

Dla średnicy $d=0{,}073\,\mathrm{m}$:

$$
 A=\pi(0{,}0365)^2\approx 0{,}00419\,\mathrm{m^2}.
$$

### 4.2. Rozłożenie siły oporu na osie

Wektor prędkości ma składowe:

$$
\vec v=(v_x,v_y),
$$

a jego wartość:

$$
 v=\sqrt{v_x^2+v_y^2}.
$$

Składowe siły oporu:

$$
 F_{Dx}=-\frac12\rho C_D A v v_x,
$$

$$
 F_{Dy}=-\frac12\rho C_D A v v_y.
$$

Z drugiej zasady Newtona:

$$
 m a_x = F_{Dx},
$$

$$
 m a_y = -mg + F_{Dy}.
$$

Ostatecznie:

$$
\boxed{
\frac{dv_x}{dt}=-\frac{\rho C_D A}{2m}\sqrt{v_x^2+v_y^2}\,v_x
}
$$

oraz:

$$
\boxed{
\frac{dv_y}{dt}=-g-\frac{\rho C_D A}{2m}\sqrt{v_x^2+v_y^2}\,v_y.
}
$$

Równania położenia pozostają następujące:

$$
\frac{dx}{dt}=v_x,\qquad \frac{dy}{dt}=v_y.
$$

### 4.3. Układ równań pierwszego rzędu

Cały model można zapisać jako układ czterech równań:

$$
\boxed{
\begin{aligned}
\dot{x} &= v_x,\\
\dot{y} &= v_y,\\
\dot{v}_x &= -k\sqrt{v_x^2+v_y^2}\,v_x,\\
\dot{v}_y &= -g-k\sqrt{v_x^2+v_y^2}\,v_y,
\end{aligned}
}
$$

gdzie:

$$
 k=\frac{\rho C_D A}{2m}.
$$

Parametr $k$ ma jednostkę $\mathrm{m^{-1}}$. Dla przyjętych
parametrów:

$$
 k=\frac{1{,}225\cdot 0{,}47\cdot 0{,}00419}{2\cdot 0{,}145}.
$$

Taki układ nie ma w ogólności prostego rozwiązania analitycznego. W
raporcie rozwiązujemy go numerycznie metodą Rungego--Kutty za pomocą
funkcji `solve_ivp` z biblioteki SciPy.

### 4.4. Warunki początkowe i warunek zakończenia obliczeń

Warunki początkowe są takie same jak w modelu idealnym:

$$
 x(0)=0,\quad y(0)=0,
$$

$$
 v_x(0)=v_0\cos\alpha,\quad v_y(0)=v_0\sin\alpha.
$$

Symulację kończymy w chwili ponownego przecięcia poziomu gruntu:

$$
 y(t)=0,\qquad t>0.
$$

W kodzie warunek ten jest realizowany jako zdarzenie numeryczne, które
zatrzymuje całkowanie przy opadającym ciele.

### 4.5. Energia mechaniczna

Energia mechaniczna jest sumą energii kinetycznej i potencjalnej:

$$
 E_m = \frac{1}{2}mv^2 + mgy.
$$

W modelu bez oporu energia mechaniczna jest zachowana. W modelu z oporem
jej pochodna wynosi:

$$
\frac{dE_m}{dt}=\vec F_D\cdot\vec v.
$$

Po podstawieniu siły oporu:

$$
\boxed{
\frac{dE_m}{dt}=-\frac12\rho C_D A v^3\le 0.
}
$$

Energia mechaniczna maleje, ponieważ część energii jest rozpraszana w
otoczeniu. Nie oznacza to naruszenia zasady zachowania energii
całkowitej --- energia jest przekazywana m.in. do powietrza w postaci
energii wewnętrznej i ruchu turbulentnego.

------------------------------------------------------------------------

## 5. Porównanie wyników i interpretacja wykresów

### 5.1. Zestawienie modeli

Po uruchomieniu skryptu wartości dla modelu z oporem są obliczane
numerycznie i wypisywane w konsoli. Model idealny daje wyniki
analityczne.

| Cecha | Bez oporu | Z oporem |
| --- | --- | --- |
| Tor ruchu | parabola | krzywa niesymetryczna |
| Składowa pozioma prędkości | stała | maleje co do wartości |
| Przyspieszenie poziome | zero | niezerowe |
| Przyspieszenie pionowe | stałe $-g$ | zależne od prędkości |
| Czas lotu | większy w tym przykładzie | krótszy lub dłuższy zależnie od parametrów |
| Zasięg | większy w tym przykładzie | mniejszy |
| Energia mechaniczna | stała | maleje |

Wartości liczbowe dla modelu oporu zależą od założonych parametrów
aerodynamicznych. Nie należy traktować ich jako uniwersalnych dla
dowolnego ciała.

### 5.2. Tor ruchu

![Porównanie torów ruchu](files/rzut_ukosny/01_tor_ruchu.png)

*Rysunek 4. Tor bez oporu jest parabolą. Opór zmniejsza składowe
prędkości i zmienia kształt toru.*

W modelu idealnym ciało wznosi się i opada symetrycznie względem
wierzchołka paraboli. Po uwzględnieniu oporu symetria zanika. Wznoszenie
i opadanie nie są lustrzanymi częściami ruchu, ponieważ prędkość oraz
siła oporu zmieniają się w czasie.

### 5.3. Prędkość i szybkość

![Wykres szybkości w funkcji czasu](files/rzut_ukosny/04_szybkosc.png)

*Rysunek 5. Szybkość jest wartością wektora prędkości, a nie jego
składową pionową.*

Bez oporu pozioma składowa prędkości jest stała. Z oporem obie składowe
zmieniają się. Warto pamiętać, że w najwyższym punkcie toru zerowa jest
tylko składowa pionowa prędkości, a nie cała prędkość.

### 5.4. Przyspieszenie

![Składowe przyspieszenia w modelu z oporem](files/rzut_ukosny/05_przyspieszenie.png)

*Rysunek 6. W modelu z oporem przyspieszenie poziome nie jest zerowe, a
pionowe nie musi być równe $-g$.*

W czasie wznoszenia opór działa w dół, zgodnie z kierunkiem przeciwnym
do prędkości. W czasie opadania opór ma składową skierowaną ku górze,
przeciwną do ruchu pionowego. Z tego powodu pionowe przyspieszenie może
mieć inną wartość niż $-g$.

### 5.5. Energia mechaniczna

![Energia mechaniczna w obu modelach](files/rzut_ukosny/06_energia.png)

*Rysunek 7. W modelu idealnym energia mechaniczna pozostaje stała,
natomiast w modelu z oporem maleje.*

Spadek energii mechanicznej wynika z pracy siły oporu. Wykres nie
przedstawia całkowitej energii układu ciało--powietrze, lecz energię
mechaniczną samego poruszającego się ciała.

### 5.6. Wartość przyspieszenia

![Wartość przyspieszenia w funkcji czasu](files/rzut_ukosny/07_wartosc_przyspieszenia.png)

*Rysunek 8. Bez oporu wartość przyspieszenia jest równa $g$. W modelu
oporu wartość ta zmienia się w trakcie lotu.*

------------------------------------------------------------------------

## 6. Ograniczenia modelu

1.  Przyjęto stałą gęstość powietrza. W rzeczywistości może ona zmieniać
    się z wysokością, temperaturą i ciśnieniem.
2.  Współczynnik $C_D$ przyjęto jako stały. W praktyce zależy on od
    kształtu ciała, orientacji i warunków przepływu.
3.  Pominięto wiatr, obrót ciała, siłę Magnusa oraz ewentualne zmiany
    orientacji.
4.  Przyjęto stałe przyspieszenie grawitacyjne
    $g=9{,}81\,\mathrm{m/s^2}$.
5.  Ciało traktowane jest jako punkt materialny, a zderzenie z podłożem
    nie jest modelowane.
6.  Parametry obiektu użyte w modelu oporu są przykładowe. Zmiana masy,
    średnicy lub współczynnika oporu zmieni wyniki.

------------------------------------------------------------------------

## 7. Kod i odtwarzanie wykresów

Wykresy zapisano w katalogu `files/rzut_ukosny/`. Skrypt
`rzut_ukosny.py` odtwarza wszystkie siedem ilustracji. Do obliczeń
numerycznych wykorzystano `solve_ivp` z pakietu SciPy, a do wykresów
Matplotlib.

### 7.1. Struktura plików

Po rozpakowaniu paczki katalog powinien wyglądać następująco:

```text
raport_rzut_ukosny/
├── raport_rzut_ukosny.md
├── rzut_ukosny.py
└── files/
    └── rzut_ukosny/
        ├── 01_tor_ruchu.png
        ├── 02_polozenie_czas.png
        ├── 03_predkosci.png
        ├── 04_szybkosc.png
        ├── 05_przyspieszenie.png
        ├── 06_energia.png
        └── 07_wartosc_przyspieszenia.png
```

### 7.2. Instalacja bibliotek

W terminalu VS Code, w katalogu projektu, wykonaj:

```bash
python -m pip install numpy scipy matplotlib
```

### 7.3. Uruchomienie skryptu

```bash
python rzut_ukosny.py
```

Program utworzy katalog `files/rzut_ukosny` (jeżeli jeszcze nie
istnieje) i zapisze w nim pliki PNG. W konsoli pojawią się również
obliczone parametry lotu dla obu modeli.

### 7.4. Kod źródłowy

Pełny, samodzielny kod znajduje się w dołączonym pliku
[`rzut_ukosny.py`](rzut_ukosny.py). Poniżej przedstawiono najważniejszy
fragment odpowiedzialny za model z oporem:

```python
def equations(t, state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)

    ax = -(K_DRAG / M) * speed * vx
    ay = -G - (K_DRAG / M) * speed * vy

    return [vx, vy, ax, ay]

def hit_ground(t, state):
    return state[1]

hit_ground.terminal = True
hit_ground.direction = -1

solution = solve_ivp(
    equations,
    (0, 30),
    [0, 0, vx0, vy0],
    events=hit_ground,
    rtol=1e-9,
    atol=1e-11,
    max_step=0.01
)
```

W powyższym fragmencie `state` zawiera kolejno $(x, y, v_x, v_y)$. Funkcja
`solve_ivp` numerycznie całkuje układ równań, a zdarzenie `hit_ground`
kończy obliczenia po ponownym osiągnięciu poziomu gruntu.

------------------------------------------------------------------------

## 8. Wnioski

1.  Rzut ukośny bez oporu powietrza można rozłożyć na dwa niezależne
    ruchy: jednostajny w poziomie i jednostajnie zmienny w pionie.
2.  Tor ruchu w modelu idealnym jest parabolą, a czas wznoszenia jest
    równy połowie całkowitego czasu lotu, jeśli start i lądowanie
    odbywają się na tej samej wysokości.
3.  Dla $v_0=36\,\mathrm{m/s}$ i $\alpha=79^\circ$ duża część prędkości
    początkowej jest skierowana pionowo. Dlatego lot trwa stosunkowo długo,
    a zasięg poziomy jest ograniczony.
4.  Uwzględnienie oporu powietrza powoduje zmianę obu składowych
    prędkości i przyspieszenia. Tor przestaje być parabolą.
5.  Siła oporu wykonuje ujemną pracę względem ruchu, dlatego energia
    mechaniczna ciała maleje.
6.  Wyniki modelu z oporem zależą od parametrów aerodynamicznych. Aby
    uzyskać wiarygodne wyniki dla konkretnego obiektu, należy znać jego
    masę, rozmiar i właściwy współczynnik oporu.
7.  Symulacja numeryczna pozwala analizować ruch także wtedy, gdy
    równania nie mają wygodnego rozwiązania analitycznego.

------------------------------------------------------------------------

## 9. Źródła i dalsza lektura

-   [OpenStax --- University Physics, Volume 1: Motion in Two Dimensions](https://openstax.org/books/university-physics-volume-1/pages/4-introduction)
-   [OpenStax --- Newton's Second Law](https://openstax.org/books/university-physics-volume-1/pages/5-3-newtons-second-law)
-   [SciPy --- `solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
-   [Matplotlib --- dokumentacja](https://matplotlib.org/stable/)

*Raport ma charakter edukacyjny. Wartości obliczeniowe są wynikami
przyjętego modelu matematycznego, a nie pomiarami laboratoryjnymi.*

Generowanie wykresów do raportu „Rzut ukośny”.
Wymagania: numpy, scipy, matplotlib.
Uruchom z katalogu zawierającego ten plik:
    python -m pip install numpy scipy matplotlib
    python rzut_ukosny.py
Wyniki trafiają do files/rzut_ukosny/.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parametry zadania
V0 = 36.0                 # m/s
ALPHA_DEG = 79.0          # stopnie
G = 9.81                  # m/s^2

# Model oporu: kula o parametrach zbliżonych do piłki tenisowej/gumowej.
M = 0.145                 # kg
RHO = 1.225               # kg/m^3, powietrze w pobliżu poziomu morza
CD = 0.47                 # współczynnik oporu (wartość modelowa)
DIAMETER = 0.073          # m
AREA = np.pi * (DIAMETER / 2)**2
K_DRAG = 0.5 * RHO * CD * AREA  # N / (m/s)^2

OUT = Path("files/rzut_ukosny")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25,
    "figure.figsize": (9, 5.4), "savefig.dpi": 180
})
BLUE, ORANGE, GREEN = "#2463A6", "#D47722", "#168578"

alpha = np.deg2rad(ALPHA_DEG)
vx0 = V0 * np.cos(alpha)
vy0 = V0 * np.sin(alpha)

# Rozwiązanie analityczne bez oporu
t_flight_ideal = 2 * vy0 / G
t_apex_ideal = vy0 / G
range_ideal = vx0 * t_flight_ideal
height_ideal = vy0**2 / (2 * G)
t_i = np.linspace(0, t_flight_ideal, 700)
x_i = vx0 * t_i
y_i = vy0 * t_i - 0.5 * G * t_i**2
vx_i = np.full_like(t_i, vx0)
vy_i = vy0 - G * t_i
speed_i = np.hypot(vx_i, vy_i)
energy_i = 0.5 * M * speed_i**2 + M * G * y_i

# Równania ruchu z oporem kwadratowym:
# F_op = -1/2 rho Cd A |v| v
def equations(t, state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)
    ax = -(K_DRAG / M) * speed * vx
    ay = -G - (K_DRAG / M) * speed * vy
    return [vx, vy, ax, ay]

def hit_ground(t, state):
    return state[1]
hit_ground.terminal = True
hit_ground.direction = -1

solution = solve_ivp(
    equations, (0, 30), [0, 0, vx0, vy0],
    events=hit_ground, rtol=1e-9, atol=1e-11, max_step=0.01
)
t_d = solution.t
x_d, y_d, vx_d, vy_d = solution.y
valid = y_d >= -1e-8
t_d, x_d, y_d, vx_d, vy_d = [a[valid] for a in (t_d, x_d, y_d, vx_d, vy_d)]
speed_d = np.hypot(vx_d, vy_d)
ax_d = -(K_DRAG / M) * speed_d * vx_d
ay_d = -G - (K_DRAG / M) * speed_d * vy_d
energy_d = 0.5 * M * speed_d**2 + M * G * y_d

def save(name):
    plt.tight_layout()
    plt.savefig(OUT / name, bbox_inches="tight")
    plt.close()

# 1. Trajektoria
plt.figure()
plt.plot(x_i, y_i, label="Bez oporu powietrza", lw=2.4, color=BLUE)
plt.plot(x_d, y_d, label="Z oporem kwadratowym", lw=2.4, color=ORANGE)
plt.scatter([range_ideal, x_d[-1]], [0, 0], s=28, color=[BLUE, ORANGE], zorder=3)
plt.xlabel("Położenie poziome x [m]"); plt.ylabel("Wysokość y [m]")
plt.title("Tor ruchu — porównanie modeli"); plt.legend(); save("01_tor_ruchu.png")

# 2. Położenie w czasie
plt.figure()
plt.plot(t_i, x_i, color=BLUE, label="x — bez oporu")
plt.plot(t_d, x_d, color=ORANGE, ls="--", label="x — z oporem")
plt.plot(t_i, y_i, color=GREEN, label="y — bez oporu")
plt.plot(t_d, y_d, color="#9467bd", ls="--", label="y — z oporem")
plt.xlabel("Czas t [s]"); plt.ylabel("Położenie [m]")
plt.title("Położenie w funkcji czasu"); plt.legend(ncol=2); save("02_polozenie_czas.png")

# 3. Składowe prędkości
plt.figure()
plt.plot(t_i, vx_i, color=BLUE, label="v_x — bez oporu")
plt.plot(t_d, vx_d, color=ORANGE, ls="--", label="v_x — z oporem")
plt.plot(t_i, vy_i, color=GREEN, label="v_y — bez oporu")
plt.plot(t_d, vy_d, color="#9467bd", ls="--", label="v_y — z oporem")
plt.axhline(0, color="black", lw=.8)
plt.xlabel("Czas t [s]"); plt.ylabel("Składowa prędkości [m/s]")
plt.title("Składowe prędkości w funkcji czasu"); plt.legend(ncol=2); save("03_predkosci.png")

# 4. Szybkość
plt.figure()
plt.plot(t_i, speed_i, color=BLUE, label="Bez oporu")
plt.plot(t_d, speed_d, color=ORANGE, label="Z oporem")
plt.xlabel("Czas t [s]"); plt.ylabel("Szybkość |v| [m/s]")
plt.title("Wartość prędkości w funkcji czasu"); plt.legend(); save("04_szybkosc.png")

# 5. Przyspieszenie w modelu z oporem
plt.figure()
plt.plot(t_d, ax_d, color=BLUE, label="a_x — z oporem")
plt.plot(t_d, ay_d, color=ORANGE, label="a_y — z oporem")
plt.axhline(-G, color="black", ls=":", label="−g (bez oporu)")
plt.axhline(0, color="gray", lw=.8)
plt.xlabel("Czas t [s]"); plt.ylabel("Przyspieszenie [m/s²]")
plt.title("Składowe przyspieszenia w modelu z oporem")
plt.legend(); save("05_przyspieszenie.png")

# 6. Energia mechaniczna
plt.figure()
plt.plot(t_i, energy_i, color=BLUE, label="Bez oporu (energia stała)")
plt.plot(t_d, energy_d, color=ORANGE, label="Z oporem (energia maleje)")
plt.xlabel("Czas t [s]"); plt.ylabel("Energia mechaniczna [J]")
plt.title("Zmiana energii mechanicznej"); plt.legend(); save("06_energia.png")

# 7. Wartość przyspieszenia
plt.figure()
plt.plot(t_d, np.hypot(ax_d, ay_d), color=GREEN, label="|a| — z oporem")
plt.plot(t_i, np.full_like(t_i, G), color=BLUE, ls="--", label="g — bez oporu")
plt.xlabel("Czas t [s]"); plt.ylabel("Wartość przyspieszenia [m/s²]")
plt.title("Wartość przyspieszenia"); plt.legend(); save("07_wartosc_przyspieszenia.png")

print("Wykresy zapisano w:", OUT.resolve())
print(f"Bez oporu: czas lotu={t_flight_ideal:.3f} s, zasięg={range_ideal:.3f} m, wysokość={height_ideal:.3f} m")
print(f"Z oporem: czas lotu={t_d[-1]:.3f} s, zasięg={x_d[-1]:.3f} m, wysokość={np.max(y_d):.3f} m")
