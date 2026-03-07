// Science commentary figure: conceptual schematic
// Compile: typst compile fig_commentary_schematic.typ fig_commentary_schematic.pdf
#set page(width: 17.8cm, height: 9.5cm, margin: (x: 0.4cm, y: 0.4cm))
#set text(font: "New Computer Modern", size: 8.5pt)

#let blue   = rgb("#1B6B93")
#let orange = rgb("#D68910")
#let red    = rgb("#C73E3A")
#let green  = rgb("#2E7D32")
#let gray   = luma(100)
#let lblue  = rgb("#E3F0F5")
#let lorange = rgb("#FDF3E4")
#let lred   = rgb("#FDEAEA")
#let lgray  = luma(245)

// Rounded panel box
#let panel(title, body, accent: blue) = {
  box(width: 100%, stroke: accent + 1pt, radius: 5pt, clip: true)[
    #box(width: 100%, fill: accent.lighten(85%), inset: (x: 8pt, y: 5pt))[
      #text(weight: "bold", size: 9pt, fill: accent)[#title]
    ]
    #box(width: 100%, inset: (x: 8pt, y: 6pt))[#body]
  ]
}

// Component node
#let node(sym, label, fill-c, stroke-c) = {
  box(inset: (x: 10pt, y: 6pt), radius: 20pt, fill: fill-c, stroke: stroke-c + 1.2pt)[
    #align(center)[
      #text(fill: stroke-c, weight: "bold", size: 11pt)[#sym] \
      #text(size: 7.5pt)[#label]
    ]
  ]
}

// Large connecting arrow
#let bigarrow = {
  align(center + horizon)[
    #text(size: 22pt, fill: luma(170), weight: "bold")[$arrow.r$]
  ]
}

// === MAIN FIGURE ===
#grid(
  columns: (5.3cm, 0.7cm, 5.3cm, 0.7cm, 5.3cm),
  align: (center + top, center + horizon, center + top, center + horizon, center + top),
  row-gutter: 0pt,

  // ── Panel 1: True DGP ──
  panel("Reality: two-component DGP", accent: blue)[
    #align(center)[
      #v(0.2em)
      #grid(
        columns: (auto, auto, auto),
        column-gutter: 0.4em,
        align: center + horizon,
        node($theta_i$, "intrinsic\nfrailty", lblue, blue),
        [
          #text(size: 13pt)[$stretch(arrow.l.r, size: #1.4em)$] \
          #text(size: 8pt, fill: gray)[$rho$]
        ],
        node($gamma_i$, "extrinsic\nsusceptibility", lorange, orange),
      )
      #v(0.5em)
      #text(size: 8pt)[
        Total hazard: \
        $mu_i (t) = underbrace(a e^(b t + theta_i), "intrinsic") + underbrace(gamma_i dot.c m_"ex", "extrinsic")$
      ]
      #v(0.4em)
      #box(inset: 4pt, radius: 3pt, fill: lgray)[
        #text(size: 7pt, fill: gray)[
          Twin concordance in $r_"MZ"$ \
          reflects _both_ genetic components
        ]
      ]
      #v(0.1em)
    ]
  ],

  bigarrow,

  // ── Panel 2: Calibration ──
  panel("Calibration: one-component model", accent: red)[
    #align(center)[
      #v(0.2em)
      #node($hat(sigma)_theta$, "sole free\nparameter", lred, red)
      #v(0.4em)
      #text(size: 8pt)[
        Fitted model has _no_ $gamma$ term. \
        All familial concordance absorbed:
      ]
      #v(0.3em)
      #box(inset: 5pt, radius: 3pt, fill: lred.lighten(50%), stroke: red.lighten(70%) + 0.5pt)[
        $hat(sigma)_theta^2 approx sigma_theta^2 + sigma_gamma^2 + 2 rho sigma_theta sigma_gamma$
      ]
      #v(0.4em)
      #box(inset: 4pt, radius: 3pt, fill: lgray)[
        #text(size: 7pt, fill: gray)[
          Calibrated $hat(sigma)_theta$ inflated by \
          +22% at default parameters
        ]
      ]
      #v(0.1em)
    ]
  ],

  bigarrow,

  // ── Panel 3: Extrapolation ──
  panel("Extrapolation: set " + $m_"ex" = 0$, accent: green.darken(10%))[
    #align(center)[
      #v(0.2em)
      #text(size: 8pt)[
        Inflated $hat(sigma)_theta$ persists after \
        removing extrinsic mortality:
      ]
      #v(0.4em)
      // Mini bar comparison
      #align(left)[
        #grid(
          columns: (2.6cm, 1.8cm),
          row-gutter: 4pt,
          align: (left, left),
          [#text(size: 7pt)[True intrinsic $h^2$:]],
          box(width: 77%, height: 11pt, fill: blue, radius: 2pt)[
            #align(right)[#text(size: 6.5pt, fill: white, weight: "bold")[#h(2pt) 0.51 #h(3pt)]]
          ],
          [#text(size: 7pt)[Estimated $hat(h)^2$:]],
          box(width: 100%, height: 11pt, fill: red, radius: 2pt)[
            #align(right)[#text(size: 6.5pt, fill: white, weight: "bold")[#h(2pt) 0.59 #h(3pt)]]
          ],
        )
      ]
      #v(0.3em)
      #text(size: 8.5pt, weight: "bold", fill: red)[
        Bias: +7.6 pp upward
      ]
      #v(0.2em)
      #text(size: 7.5pt)[
        Persists under 5-target calibration (+6.4 pp) \
        and with cutoff = 0 (+3.6 pp)
      ]
      #v(0.1em)
    ]
  ],
)

// ── Bottom annotation ──
#v(0.3em)
#align(center)[
  #box(inset: (x: 12pt, y: 5pt), radius: 4pt, fill: lgray, stroke: luma(210) + 0.5pt)[
    #grid(
      columns: (auto, 2em, auto, 2em, auto),
      align: center + horizon,
      column-gutter: 0pt,
      text(size: 8pt, fill: orange, weight: "bold")[$rho > 0$: bias upward #sym.arrow.t],
      [],
      text(size: 8pt, fill: gray)[#sym.bar.v],
      [],
      text(size: 8pt, fill: blue, weight: "bold")[$rho < 0$: bias reverses #sym.arrow.b],
    )
    #v(1pt)
    #text(size: 7pt, fill: gray)[Sign reversal under negative pleiotropy — a falsifiable prediction confirmed by simulation]
  ]
]
