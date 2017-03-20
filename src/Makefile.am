bin_PROGRAMS = lavax

# -------------------- JAVA --------------------
JAVAROOT = .

dist_noinst_JAVA = lavax.java Vec3.java Particle.java
# dist_noinst_SOURCES = lavax.conf predictor.in
jar_SOURCES = lavax.conf predictor.in
AM_JAVACFLAGS = -target 1.7 -source 1.7

jarfile = lavax.jar
jardir = $(datadir)/lavax

$(jarfile): $(dist_noinst_JAVA) $(jar_SOURCES)
	jar cfev $@ lavax *.class $(jar_SOURCES)

jar_DATA = $(jarfile)

lavax: $(dist_noinst_JAVA)
	echo "java -jar "$(jardir)"/lavax.jar" > lavax
	chmod +x lavax

CLEANFILES = $(jarfile) *.class
