library (dplyr)
library (cripto2)
library(readxl)

# Preluați o listă de monede criptomonede:
  
monede <- crypto_list(only_active = TRUE)
#Funcția crypto_list din pachetul crypto2 este utilizată pentru a prelua o listă de monede criptomonede.
#Argumentul only_active = TRUE specifică faptul că numai criptomonedele active ar trebui incluse în listă.
#Încărca o listă cu Top 100 de monede pe 2023-09-09 dintr-un fișier Excel:
  
top100 <- read_xlsx("top100_2023-09-09.xlsx")
# Funcția read_xlsx din pachetul readxl este utilizată pentru a citi un fișier Excel numit „top100_2023-09-09.xlsx” și pentru a stoca conținutul acestuia în cadrul de date top100.
# Acest fișier Excel conține probabil o listă cu primele 100 de monede criptomonede la data specificată.
# Alătură informațiile despre monede cu Lista Top 100:
  
top100_coins <- top100 %>%
left_join(monede, c("Nume" = "nume", "simbol"))

# Operatorul %>% (conducta) este folosit pentru a înlănțui mai multe operații într-o conductă de date. În acest caz, permite efectuarea pașii de manipulare a datelor secvenţial.
# Funcția left_join din pachetul dplyr este utilizată pentru a combina  cadrul de date top100 cu cadrul de date monede pe baza coloanelor care se potrivesc.
# Coloanele utilizate pentru alăturare sunt specificate ca „Nume” = „nume” și „simbol”.
# Această operațiune adaugă informații despre monede din cadrul de date pentru monede în cadrul de date de top 100 pe baza numelor și simbolurilor de monede care se potrivesc.
# După rularea acestui cod, cadrul de date top100_coins va conține un set de date îmbinat cu informații despre primele 100 de criptomonede la data specificată, inclusiv detalii din cadrul de date monede (de exemplu, nume, simbol). 

# Deci este un cod R conceput pentru a descărca datele serii cronologice pentru primele 100 de monede criptomonede dintr-un interval de date specificat și pentru a salva datele rezultate într-un fișier RData. Iată o explicație pas cu pas a codului:
  
# Inițializare   variabilă de buclă i:
  
i = 1
# Această linie inițializează o variabilă de contor de buclă i la 1. Această variabilă va fi folosită pentru a itera prin fiecare rând din cadrul de date top100_coins.
# For Loop pentru a descărca datele din seria temporală:
  
for (i in 1:nrow(top100_coins)) {
# Această buclă for iterează prin rândurile cadrului de date top100_coins. nrow(top100_coins) returnează numărul de rânduri din cadrul de date, care determină numărul de iterații.
# Extrage simbolul monedei și afișați un mesaj:
    simbol_ <- top100_coins$simbol[i]
  mesaj(simbol_)}

# În interiorul buclei, extrage simbolul monedei criptomonede din cadrul de date top100_coins și îl stochează în variabila simbol_.
# Funcția de mesaj este utilizată pentru a afișa valoarea simbol_ în consola R, ceea ce poate fi util pentru urmărirea progresului.
# Descărcarea datelor seriei temporale pentru monedă:
      
quot_ <- crypto_history(monede[i,],
                        start_date = "20130501",
                        end_date = "20230910")
# Acest segment de cod folosește funcția crypto_history pentru a descărca datele seriilor temporale pentru moneda criptomonedă curentă specificată de coins[i,].
# Specifică o dată de începere „2013 05 01” și o dată de încheiere „2023 09 10” pentru datele seriei cronologice.
# Adăugați sau inițializați citations_all Cadrul de date:
      
if (!exists("quotations_toate")) citations_all <- quot_ else
  quotations_all <- rbind(quotations_all, quot_)
# Acest segment de cod verifică dacă cadrul de date quotations_all există. Dacă nu există, îl inițializează cu cadrul de date quot_. Dacă există, adaugă datele quot_ la cadrul de date existent quotations_all folosind rbind.
# Acest lucru are ca rezultat un cadru de date cumulativ care conține date de serie cronologică pentru toate monedele criptomonede procesate până acum.
# Eliminare variabilele temporare:
      
rm(simbol_, quot_)

    # După descărcarea și adăugarea datelor din seria temporală pentru o monedă, elimină variabilele temporare simbol_ și quot_ din mediul R pentru a elibera memorie.
# Salvare date finale într-un fișier RData:
      
save(list = "quotations_all", file = "quotations_top100.RData")
# În cele din urmă, după procesarea tuturor monedelor, salvează cadrul de date quotations_all într-un fișier RData numit „quotations_top100.RData”. Acest fișier va conține datele combinate ale seriei cronologice pentru toate primele 100 de monede criptomonede.
# În rezumat, acest cod iterează prin primele 100 de monede criptomonede, descarcă datele serii cronologice pentru fiecare monedă, anexează datele la un cadru de date cumulate și apoi salvează datele combinate într-un fișier RData pentru analize ulterioare.
  
    
    