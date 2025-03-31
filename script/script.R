#Tarea 3 

#cargar paquetes 
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)
library(dplyr)

#cargar objeto 
data("dietswap", package = "microbiome")
ps <- dietswap
ps


#- ¿Cuántas muestras y taxones contiene el objeto?
#tiene 222 muestras y 130 taxones
#- ¿Qué variables están disponibles en los metadatos de las muestras?
#poner en la consola View(sample_data(ps))
#sex: masuculino o femenino, nationality: AMM y AFR, group (dieta):HE, DI y ED, time point :1,2,3,4,5,6, time point within group: 1,2 subject, BMI (indice de masa corporal): lean, overweight, obese, sample


#crea la curva de rarefraccion.
otu_table(ps)
otu_table(ps)-> ab
ab <- as.data.frame(t(ab)) #intercambiamos filas y columnas
View(ab)

### Curvas de rarefacción
rarecurve(ab,col = rainbow (length(ab)),  step  =  11000, cex = 0.5)

pdf("figuras/curva_rarefaccion.pdf", height = 11, width =15)
rarecurve(ab, col = rainbow(nrow(ab)), step = 11000, cex = 0.5)
dev.off()



#1.  ¿Qué indican estas curvas? indican si el muestreo es suficiente para representar la diverdidad de la muestra
#2.  ¿Hay muestras que deberían descartarse por bajo conteo?
#hay una muestea en partixular que tienen un tamaño de muestra menor en comparacion con todas las demas muestras en este caso la muestra 56 por lo que podria considerar eliminarla. 


#####################
#diversidad alfa 
pdf("figuras/diversidad_alfa.pdf", height = 11, width =15)
plot_richness(ps,color="group", x="group" ,title = "diversidad alfa", measures =c("Observed", "Shannon", "Simpson")) + geom_boxplot(aes(fill = group), alpha=.7) + scale_color_manual(values = c("#FFE1FF", "#b2df8a", "#cae1ff")) + scale_fill_manual(values = c("#FFE1FF", "#b2df8a", "#cae1ff"))
dev.off()

pdf("figuras/diversidad_alfa2.pdf", height = 11, width =15)
plot_richness(ps,color="nationality", x="nationality" ,title = "diversidad alfa2", measures =c("Observed", "Shannon", "Simpson")) + geom_boxplot(aes(fill = nationality), alpha=.7) + scale_color_manual(values = c("#FFE1FF", "#b2df8a")) + scale_fill_manual(values = c("#FFE1FF", "#b2df8a"))
dev.off()

#HE home environment entorno habitual dietas normales antes
#ED dia de endoscopia 
#DI dietary intervention dieta especifica 
#para la diversidad en relacion a la nacionalidad los african american tienen una riqueza promedio menor a los africanos rurales, en cuando a shannon indica una uniformidad mayor en las muestras de african american y en cuento a simpson existe una dominancia en las muestras provenientes de african american teniendo un promedio superio a 0.8 en comabio las muestras de arfricanos rurales muestran menor dominancia con promedio menor a .8.
#para la diversidad en cuanto a la dieta, riqueza existe una riqueza promedio mayor en la dieta normal, en shanon  indican una abundancia poco equitativa aunque aumenta un poco durante la la intervencion, deacuerdo a simpson existe una dominancia de pocas especies por lo que hay baja diversidad sin embargo la dominancia disminuye un poco en la intervencion dietetica.
#existe una diferencia notoria entre las muestras segun la nacionalidad esto puede observarse especialmente en indice de simpson donde la dominancia de pocas especies es notablemente mayor.
#en cuento a la dieta si bien se pueden observar cambios antes y despues de la intervencion estos no son tan evidentes como lo es la nacionalidad. 

### Filtrado y transformación

#Aplica un filtrado para quedarte solo con los géneros más abundantes (por ejemplo, los que tienen más del 0.1% de abundancia relativa en al menos 10% de las muestras).
filtrado_ps<- ps %>% transform_sample_counts(function (x)x/sum(x)) %>% # se transforma para tener abundancias relativas al dividir cada valor entre el total de cada muestra
  filter_taxa(function(x) mean(x>0.001) >0.1, prune = TRUE) #se filtran x>0.001 abundancia relativa mayor a 0.001 taxon, mean(x>0.001) muestras donde sea superior a 0.1 en todas las muestras
View(filtrado_ps)

otu_filt<- otu_table (filtrado_ps)
View(otu_filt)

sample_data(otu_filt)
### Diversidad beta

#Realiza una ordención PCoA utilizando distancia Bray-Curtis. 
#Usa `ordinate()` y `plot_ordination()`.

div_bet<- ordinate ( filtrado_ps, method = "PCoA", distance = "bray")
View(div_bet)

pdf("figuras/ordenacion.pdf", height = 11, width =15)
plot_ordination(filtrado_ps, div_bet, color = "nationality" )
dev.off()
#Responde:

#- ¿Los grupos se separan visiblemente? sí los grupos se separan notablemente deacuerdo a su nacionalidad   quedando las muestras procedientes de personas afroamericanas en la aprte superior derecha e inferior izquierda, en cambio las muestras de africa rural se ubican en la parte inferior derecha y un como arriba derecha tambien en menor proporcion abajo izquierda. 
# - ¿Qué podría estar causando esas diferencias? la abundancia y riqueza de especies entre las personas afroamericanas y africanas rurales son distintas debido a la dieta y entorno en el que se desarrollan las personas


### Gráfica *Rank-Abundance*

#Crea una gráfica de abundancia-rango para observar la dominancia de taxones.


library(dplyr)

filtrado_ps #abundancia relativa

obj_f<- psmelt(filtrado_ps) #hacer un obj

ab_tx<- aggregate(Abundance ~ OTU + Phylum, data = obj_f, mean) #calcualr media para cada tax

#grafica 
pdf("figuras/Rank_abundance_puntos.pdf", height = 11, width =15)
ggplot(ab_tx, aes(x = reorder(OTU, -Abundance), y = Abundance)) +
  geom_point(aes(color = Phylum), size = 3) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10() +
  labs(
    title = "Rank abundance puntos", 
    x = "OTU (Ordenados por Abundancia)",
    y = "abundancia",
    fill = "Phylum")
dev.off()

pdf("figuras/Rank_abundance_barras.pdf", height = 11, width =15)
ggplot(ab_tx, aes(x = reorder(OTU, -Abundance), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", alpha = 0.8) + 
  labs(
    title = "rank abundance barras",
    x = "OTU (Ordenados por Abundancia)",
    y = "abundancea",
    fill = "Phylum"
  ) 
dev.off()

#Interpreta:

#- ¿Qué tan dominada está la comunidad por pocos taxones y cuáles son? son 2 phylum los dominantes bacterioidetes y firmicutes especialmente bacteroidetes tiene una abundancia mucho mayor
#- ¿Qué tipo de distribución parece seguir? una distribucion zipf

############
### Gráficas apiladas de abundancia por taxón
library(phyloseq)
library(microbiome)
library(seriation)

part_ps<- filtrado_ps %>% 
  tax_glom("Phylum")  #agurpamos por phylum

##
part_ps <- tax_glom(filtrado_ps, "Phylum")
da_m<- psmelt(part_ps)

pdf("figuras/Gráficas_apiladas_de_abundancia.pdf", height = 11, width =15)
ggplot(da_m, aes(x=Sample, y=Abundance, fill= Phylum))+  
  geom_bar(stat = "identity", position = "stack")+
  labs(title ="grafica apilada de abundancia", x="muestra", y="abundancia", fill= "Phylum" )
dev.off()

#- ¿Hay algún phylum que domine? los phylums dominantes es bacteroidetes y firmicutes
# - ¿Se observan diferencias entre grupos?  en general todos los grupos tienen una dominancia de bacteroidetes y firmicutes y poca precencias de otros phylums como actinobacterias y proteobacterias

### Exportar resultados

#Exporta alguna tabla de abundancias o métricas de diversidad a CSV.

diversidad_a<- estimate_richness(ps, measures = c("observed", "Shannon", "Simpson"))
write.csv(diversidad_a, file = "csv/diversidad_alfa.csv", row.names = TRUE)

###  `GlobalPatterns`
#Usaremos el dataset `GlobalPatterns` incluido en `phyloseq`. 
#Contiene 26 muestras de diversos ambientes.

data("GlobalPatterns")
gp <- GlobalPatterns
#filtrar
gp_f<- gp %>% 
  filter_taxa(function(x) sum (x>5) > 0.2 * length(x), prune =TRUE)#almenos 5 lecturas en 20%, aquellas que no cumplen se eliminan

#aglomerar 
ag<- tax_glom(gp_f, "Family")

abr2<- transform_sample_counts(ag, function (x) x/ sum (x)*100)

#subset
f_sfs<- subset_samples( abr2, SampleType %in% c("Soil", "Feces", "Skin"))
print(f_sfs)
#sin abun relativa 
f_sfs_sar<- subset_samples( ag, SampleType %in% c("Soil", "Feces", "Skin"))
#dim
print(f_sfs_sar)

#- Calcular 3 índices de diversidad alfa (`Shannon`, `Simpson`, `Observed`)

diversidad_a<- estimate_richness(f_sfs_sar, measures = c("Observed", "Shannon", "Simpson"))
print(diversidad_a) #no quiere abundancias relativas
#- Crear boxplots comparativos de los índices entre tipos de muestra

library(ggplot2)
diversidad_a$SampleType <- sample_data(f_sfs_sar)$SampleType #agregar sample type para relacionar con soli, feces y skin

# Observed
pdf("figuras/Box plot D. observada.pdf", height = 11, width =15)
ggplot(diversidad_a, aes(x = SampleType, y = Observed, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Box plot D. observada")
dev.off()

pdf("figuras/Box plot D.shannon.pdf", height = 11, width =15)
ggplot(diversidad_a, aes(x = SampleType, y = Shannon, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Box plot D. shannon")
dev.off()


pdf("figuras/Box plot D.simpson .pdf", height = 11, width =15)
ggplot(diversidad_a, aes(x = SampleType, y = Simpson, fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Box plot D. simpson")
dev.off()

#- Realizar prueba estadística (Kruskal-Wallis) para diferencias entre grupos

k_o <- kruskal.test(Observed ~ SampleType, data = diversidad_a)
print(k_o)


k_si<- kruskal.test(Simpson ~ SampleType, data = diversidad_a)
print(k_si)

k_sha <- kruskal.test(Shannon ~ SampleType, data = diversidad_a)
print(k_sha)


### Curvas de Rango-Abundancia 

#Crear gráficas de rango-abundancia para cada tipo de muestra
#Usar escala log10 en Y
#Comparar patrones entre ambientes

obj_f<- psmelt(f_sfs_sar) #convertir de phy a data

ab_tx_muestras <- aggregate(Abundance ~ OTU + Phylum ,data = obj_f, mean)

ggplot(ab_tx_muestras,aes(x=reorder(OTU,-Abundance),y=Abundance, fill = Phylum, color = Phylum)) +
  geom_point()+scale_y_log10()

### Perfil taxonómico 

#- Crear gráfico apilado de abundancia a nivel de Phylum
#- Mostrar solo los 5 phyla más abundantes
#- Agrupar por tipo de muestra
#- Usar `facet_wrap` para comparar ambientes

d_b<- ordinate (f_sfs_sar, method = "PCoA", distance = "bray")
View(d_b)
class(d_b)

#grafica

plot_ordination(f_sfs_sar, d_b, color = "SampleType") +
  stat_ellipse(aes(color = SampleType), level = 0.95)
dev.off()


library(vegan)
bray_dist <- phyloseq::distance(f_sfs_sar, method = "bray")
md<- as(sample_data(f_sfs_sar), "data.frame")

per<- adonis2(bray_dist ~ SampleType,data = md )
print (per)
