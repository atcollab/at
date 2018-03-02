% NAFFLIB MATLAB to NAFF library
%
%  INPUTS
%  1. Real part
%  2. Imaginary part
%  3. Window type
%  4. Number of frequencies
%  5. Debug 0 or 1
%
%  OUPUTS
%  1. Fundamental frequency vector
%  2. Amplitude vector
%  3. Phase vector
%
%  EXAMPLE
%  1. [frequency amplitude phase] = nafflib(Y, Yp, WindowType,nfreq,DebugFlag);


%% Laurent S. Nadolski, Synchrotron SOLEIL

% Calls modnaff from: 
% /* MODNAFF.C version 0.96                                                 */
% /* Astronomie et systemes dynamiques                                      */
% /* J. LASKAR Bureau des Longitudes, Paris (version originale fortran)     */
% /* M. GASTINEAU, portage en C, 03/09/98                                   */
% /* v0.96 M. GASTINEAU 09/09/98 : modification dans la fonction naf_ztder  */
% /*                               dans les boucles de I en I-1.            */
% /* v0.96 M. GASTINEAU 09/11/98 : remplacement des variables PI2 par PICARRE*/
% /*                               car PI2 est un define sous LINUX.        */
% /* v0.96 M. GASTINEAU 01/12/98 : utilisation des fonctions  complexes     */
% /*                               inlines.                                 */
% /* v0.96 M. GASTINEAU 07/12/98 : modification dans naf_frefin car bug lors*/
% /*                               de l'optimisation sur Mac.               */
% /* v0.96 M. GASTINEAU 18/12/98 : modification de naf_iniwin               */
% /* v0.96 M. GASTINEAU 06/01/99 : ajout de la liberation de la liste des   */
% /*                        fenetres dans naf_cleannaf et naf_cleannaf_notab*/
% /* v0.96 M. GASTINEAU 12/01/99 : correction dans le cas ou ICPLX=0 dans   */
% /*                        naf_modfre et naf_gramsc.                       */
% /*v0.96 M. GASTINEAU 14/01/99 : ajout du support des fenetres dans        */
% /*                        naf_mfttab et naf_fftmax.                       */
% /*v0.97 M. GASTINEAU 26/05/99 : correction bug (0.97/99/05/26/A) dans     */
% /*                        naf_mftnaf                                      */
% /* v0.97 M. GASTINEAU 27/05/99 : correction bug (0.97/99/05/27/A) dans le */
% /*                  cas ou ICPLX=0 et FS=0 dans naf_modfre et naf_gramsc. */



