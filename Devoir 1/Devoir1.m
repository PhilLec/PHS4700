function [pcm MI aa] = Devoir1(pos, theta, wz, Force)
   
   # matrice de rotation (autour de y)
   ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
   
   # masses
   m_jambe = (pi * 0.06^2 * 0.75) * 1052;
   m_tronc = (pi * 0.15^2 * 0.7) * 953;
   m_cou = (pi * 0.04^2 * 0.1) * 953;
   m_tete = (4 * pi * 0.1^3 / 3) * 1056;
   m_bras = (pi * 0.03^2 * 0.75) * 1052;
   m_totale = 2 * m_jambe + m_tronc + m_cou + m_tete + 2 * m_bras;
   
   # centres de masses dans le référentiel du patineur sans rotation
   cm_jambe_1 = [-0.1; 0; 0.75/2];
   cm_jambe_2 = [0.1; 0; 0.75/2];
   cm_tronc = [0; 0; 0.75 + 0.7/2];
   cm_cou = [0; 0; 0.75 + 0.7 + 0.1/2];
   cm_tete = [0; 0; 0.75 + 0.7 + 0.1 + 0.1];
   cm_bras_1 = [-(0.15 + 0.03); 0; 0.75 + 0.7 - 0.75/2];
   cm_bras_2 = [0.15 + 0.75/2; 0; 0.75 + 0.7 + 0.03/2];
   
   # CM considérant theta
   cm_jambe_1 = ry * cm_jambe_1;
   cm_jambe_2 = ry * cm_jambe_2;
   cm_tronc = ry * cm_tronc;
   cm_cou = ry * cm_cou;
   cm_tete = ry * cm_tete;
   cm_bras_1 = ry * cm_bras_1;
   cm_bras_2 = ry * cm_bras_2;
   
   # CM en considérant pos
   cm_jambe_1 += pos;
   cm_jambe_2 += pos;
   cm_tronc += pos;
   cm_cou += pos;
   cm_tete += pos;
   cm_bras_1 += pos;
   cm_bras_2 += pos;
   
   pcm = (m_jambe * cm_jambe_1 + m_jambe * cm_jambe_2 + m_tronc * cm_tronc + m_cou * cm_cou + m_tete * cm_tete + m_bras * cm_bras_1 + m_bras * cm_bras_2) / m_totale;
   
   # moments d'inertie
   mi_jambe = [m_jambe * 0.06^2 / 4 + m_jambe * 0.75^2 / 12; m_jambe * 0.06^2 / 4 + m_jambe * 0.75^2 / 12; m_jambe * 0.06^2 / 2];
   mi_tronc = [m_tronc * 0.15^2 / 4 + m_tronc * 0.7^2 / 12; m_tronc * 0.15^2 / 4 + m_tronc * 0.7^2 / 12; m_tronc * 0.15^2 / 2];
   mi_cou = [m_cou * 0.04^2 / 4 + m_cou * 0.1^2 / 12; m_cou * 0.04^2 / 4 + m_cou * 0.1^2 / 12; m_cou * 0.04^2 / 2];
   mi_tete = [2 * m_tete * 0.1^2 / 5; 2 * m_tete * 0.1^2 / 5; 2 * m_tete * 0.1^2 / 5];
   mi_bras_para = [m_bras * 0.03^2 / 4 + m_bras * 0.75^2 / 12; m_bras * 0.03^2 / 4 + m_bras * 0.75^2 / 12; m_bras * 0.03^2 / 12];
   mi_bras_perp = [m_bras * 0.03^2 / 12; m_bras * 0.03^2 / 4 + m_bras * 0.75^2 / 12; m_bras * 0.03^2 / 4 + m_bras * 0.75^2 / 12];
   
   # MI en considérant theta
   mi_jambe = ry * mi_jambe;
   mi_tronc = ry * mi_tronc;
   mi_cou = ry * mi_cou;
   mi_tete = ry * mi_tete;
   mi_bras_para = ry * mi_bras_para;
   mi_bras_perp = ry * mi_bras_perp;
   
   # MI en considérant pos
   mi_jambe_1_g = getMMI(pos, mi_jambe, m_jambe, cm_jambe_1);
   mi_jambe_2_g = getMMI(pos, mi_jambe, m_jambe, cm_jambe_2);
   mi_tronc_g = getMMI(pos, mi_tronc, m_tronc, cm_tronc);
   mi_cou_g = getMMI(pos, mi_cou, m_cou, cm_cou);
   mi_tete_g = getMMI(pos, mi_tete, m_tete, cm_tete);
   mi_bras_1_g = getMMI(pos, mi_bras_para, m_bras, cm_bras_1);
   mi_bras_2_g = getMMI(pos, mi_bras_perp, m_bras, cm_bras_2);
   
   MI = mi_jambe_1_g + mi_jambe_2_g + mi_tronc_g + mi_cou_g + mi_tete_g + mi_bras_1_g + mi_bras_2_g;
   
   # rayon r_f
   zf = 0.75 + 0.70 + 0.10 + 0.10;
   rf = [0, 0.10, zf]';
   
   # rayon en considérant theta
   rf = ry * rf;
   
   # calcul du torque
   t = cross(rf, Force);
   
   # calcul de l'accélération angulaire
   aa = inv(MI) * t;
   
   # aa = [0; 0; 0];
   
 end