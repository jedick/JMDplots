% extracted from CHNOSZ's macros.Rd 20190929
%\newcommand{\ZC}{\ifelse{latex}{\eqn{Z_\mathrm{C}}}{\ifelse{html}{\out{<I>Z</I><sub>C</sub>}}{ZC}}}
\newcommand{\logK}{\ifelse{latex}{\eqn{\log K}}{\ifelse{html}{\out{log<i>K</i>}}{log K}}}
% 20190930
\newcommand{\nH2O}{\ifelse{latex}{\eqn{n\mathrm{H_{2}O}}}{\ifelse{html}{\out{<i>n</i>H<sub>2</sub>O}}{nH2O}}}
% 20191005
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}
% 20191020
\newcommand{\pKa}{\ifelse{latex}{\eqn{\mathrm{p}K_\mathrm{a}}}{\ifelse{html}{\out{p<i>K</i><sub>a</sub>}}{pKa}}}
% 20191021
\newcommand{\T}{\ifelse{latex}{\eqn{T}}{\ifelse{html}{\out{<I>T</I>}}{T}}}
\newcommand{\Topt}{\ifelse{latex}{\eqn{T_\mathrm{opt}}}{\ifelse{html}{\out{<I>T</I><sub>opt</sub>}}{Topt}}}
% 20191024
\newcommand{\O2}{\ifelse{latex}{\eqn{\mathrm{O_{2}}}}{\ifelse{html}{\out{O<sub>2</sub>}}{O2}}}
% 20200204
\newcommand{\nAA}{\ifelse{latex}{\eqn{n_{\mathrm{AA}}}}{\ifelse{html}{\out{<i>n</i><sub>AA</sub>}}{nAA}}}
\newcommand{\nO2}{\ifelse{latex}{\eqn{n\mathrm{O_{2}}}}{\ifelse{html}{\out{<i>n</i>O<sub>2</sub>}}{nO2}}}
% 20201015
% https://stackoverflow.com/questions/595365/how-to-render-narrow-non-breaking-spaces-in-html-for-windows
\newcommand{\logfO2}{\ifelse{latex}{\eqn{\log f\mathrm{O_{2}}}}{\ifelse{html}{\out{log&#x202F;<i>f</i>O<sub>2</sub>}}{logfO2}}}
\newcommand{\logaH2O}{\ifelse{latex}{\eqn{\log a\mathrm{H_{2}O}}}{\ifelse{html}{\out{log&#x202F;<i>a</i>H<sub>2</sub>O}}{logaH2O}}}
% 20210205
\newcommand{\degC}{\ifelse{latex}{\eqn{^{\circ}}C}{\ifelse{html}{\out{&deg;}C}{°C}}}
% 20210401
\newcommand{\Hplus}{\ifelse{latex}{\eqn{\mathrm{H^{+}}}}{\ifelse{html}{\out{H<sup>+</sup>}}{H+}}}
% 20210816
\newcommand{\CO2}{\ifelse{latex}{\eqn{\mathrm{CO_{2}}}}{\ifelse{html}{\out{CO<sub>2</sub>}}{CO2}}}
\newcommand{\H2}{\ifelse{latex}{\eqn{\mathrm{H_{2}}}}{\ifelse{html}{\out{H<sub>2</sub>}}{H2}}}
\newcommand{\NH4plus}{\ifelse{latex}{\eqn{\mathrm{NH_{4}^{+}}}}{\ifelse{html}{\out{NH<sub>4</sub><sup>+</sup>}}{NH4+}}}
\newcommand{\H2S}{\ifelse{latex}{\eqn{\mathrm{H_{2}S}}}{\ifelse{html}{\out{H<sub>2</sub>S}}{H2S}}}
\newcommand{\CH4}{\ifelse{latex}{\eqn{\mathrm{CH_{4}}}}{\ifelse{html}{\out{CH<sub>4</sub>}}{CH4}}}
% 20220420
\newcommand{\logaH2}{\ifelse{latex}{\eqn{\log a\mathrm{H_{2}}}}{\ifelse{html}{\out{log&#x202F;<i>a</i>H<sub>2</sub>}}{logaH2}}}
% Subscripts and superscripts 20220520
\newcommand{\s}{\ifelse{latex}{\eqn{_{#1}}}{\ifelse{html}{\out{<sub>#1</sub>}}{#1}}}
\newcommand{\S}{\ifelse{latex}{\eqn{^{#1}}}{\ifelse{html}{\out{<sup>#1</sup>}}{^#1}}}
% To get mu in LaTeX 20220621
\newcommand{\umol}{\ifelse{latex}{\eqn{\mu}mol}{\ifelse{html}{\out{&mu;mol}}{μmol}}}
% \ZC changed to \Zc (and canprot::ZCAA() renamed to canprot::Zc()) 20230617
\newcommand{\Zc}{\ifelse{latex}{\eqn{Z_\mathrm{C}}}{\ifelse{html}{\out{<I>Z</I><sub>C</sub>}}{Zc}}}
% 20240206
\newcommand{\Cp}{\ifelse{latex}{\eqn{C_P}}{\ifelse{html}{\out{<I>C<sub>P</sub></I>}}{Cp}}}
