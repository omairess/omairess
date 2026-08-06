'use strict';

/* BSRT translations.
 *
 * SHARED FILE — a byte-identical copy lives at BSRT-desktop/renderer/i18n.js.
 *
 * TRANSLATION STATUS: English is the source. The French, Dutch and German
 * strings were produced for this app and have NOT been through a formal
 * translation-verification process. The Karolinska Sleepiness Scale anchors in
 * particular are part of a validated instrument — if you publish, check them
 * against an officially validated translation for your language and correct
 * this file if they differ. The language actually used is exported with every
 * trial (`language`), so it is always visible which wording a participant saw.
 *
 * Applying translations: elements carry data-i18n="key" for text content and
 * data-i18n-ph="key" for placeholders. Call applyTranslations() after changing
 * language. Anything without a data-i18n attribute stays in English.
 */

(function (root, factory) {
  var api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.BSRTi18n = api;
})(typeof self !== 'undefined' ? self : this, function () {

  var LANGUAGES = [
    { code: 'en', label: 'English' },
    { code: 'fr', label: 'Français' },
    { code: 'nl', label: 'Nederlands' },
    { code: 'de', label: 'Deutsch' }
  ];

  /*
   * Karolinska Sleepiness Scale (Akerstedt & Gillberg, 1990), 9-point form.
   * Anchors are given for every step, which is the form most commonly used in
   * sleep research.
   */
  var KSS = {
    en: ['Extremely alert', 'Very alert', 'Alert', 'Rather alert',
         'Neither alert nor sleepy', 'Some signs of sleepiness',
         'Sleepy, but no effort to keep awake', 'Sleepy, some effort to keep awake',
         'Very sleepy, great effort to keep awake, fighting sleep'],
    fr: ['Extrêmement alerte', 'Très alerte', 'Alerte', 'Plutôt alerte',
         'Ni alerte ni somnolent', 'Quelques signes de somnolence',
         'Somnolent, mais sans effort pour rester éveillé',
         'Somnolent, un certain effort pour rester éveillé',
         'Très somnolent, gros effort pour rester éveillé, lutte contre le sommeil'],
    nl: ['Extreem alert', 'Zeer alert', 'Alert', 'Tamelijk alert',
         'Niet alert, niet slaperig', 'Enkele tekenen van slaperigheid',
         'Slaperig, maar geen moeite om wakker te blijven',
         'Slaperig, enige moeite om wakker te blijven',
         'Zeer slaperig, grote moeite om wakker te blijven, vecht tegen de slaap'],
    de: ['Äußerst wach', 'Sehr wach', 'Wach', 'Ziemlich wach',
         'Weder wach noch schläfrig', 'Erste Anzeichen von Schläfrigkeit',
         'Schläfrig, aber ohne Mühe wach zu bleiben',
         'Schläfrig, etwas Mühe wach zu bleiben',
         'Sehr schläfrig, große Mühe wach zu bleiben, kämpfe gegen den Schlaf']
  };

  var STRINGS = {
    en: {
      'app.subtitle': 'Behavioral Sleep Resistance Task and Psychomotor Vigilance Task.',
      'lang.label': 'Language',

      'participant.heading': 'Participant',
      'participant.id': 'Participant ID',
      'participant.name': 'Name',
      'participant.birth': 'Birth date',
      'participant.education': 'Educational level',
      'participant.address': 'Address',
      'participant.session': 'Session / condition',
      'participant.trial': 'Trial number',
      'participant.privacy': 'Name, address and birth date are directly identifying. They are stored on this device and written into every exported file. For research use, prefer a pseudonymous participant ID and keep the identity key separately under your data-protection plan. Date and time of testing are captured automatically.',

      'task.heading': 'Task',
      'task.mode': 'Mode',
      'task.mode.bsrt': 'BSRT / OSLER — fixed interval',
      'task.mode.pvt': 'PVT — variable interval',
      'task.isi': 'Stimulus interval (ms)',
      'task.isiSet': 'Intervals (seconds, comma separated)',
      'task.block': 'Block length (s)',
      'task.seed': 'Schedule seed (blank = random)',

      'protocol.heading': 'Protocol',
      'protocol.stim': 'Stimulus duration / hit window (ms)',
      'protocol.max': 'Maximum duration (min)',
      'protocol.lapse': 'Lapse threshold (ms)',
      'protocol.criterion': 'Consecutive misses = sleep onset',
      'protocol.criterionOn': 'End the test early on the sleep-onset criterion',
      'protocol.alarm': 'Sound an alarm when the sleep-onset criterion is reached',
      'protocol.fullscreen': 'Fullscreen during the task',
      'protocol.explain': 'A response within the hit window is a hit; a slow hit (beyond the lapse threshold) is a lapse. A response later than the hit window is a miss — its raw reaction time is still recorded, but it does not reset the consecutive-miss run. A response under the false-start threshold is a false start.',

      'kss.heading': 'Karolinska Sleepiness Scale',
      'kss.when': 'Administer the KSS',
      'kss.when.none': 'Not at all',
      'kss.when.before': 'Before the task',
      'kss.when.after': 'After the task',
      'kss.when.both': 'Before and after',
      'kss.question': 'How sleepy do you feel right now?',
      'kss.instruction': 'Choose the statement that best describes how you feel at this moment.',
      'kss.beforeTitle': 'Before the task',
      'kss.afterTitle': 'After the task',
      'kss.continue': 'Continue',

      'corrections.heading': 'Corrections',
      'corrections.type': 'Correction',
      'corrections.none': 'None',
      'corrections.falseStarts': 'Remove false starts only',
      'corrections.outliers': 'Remove outliers only',
      'corrections.both': 'Remove both',
      'corrections.threshold': 'False-start threshold (ms)',
      'corrections.sd': 'Outlier cut-off (SD)',

      'instructions.heading': 'Instructions to read to the participant',
      'instructions.bsrt': '“A red light will flash every few seconds. Each time you see it, press the space bar. Please try to stay awake, but do not do anything else to keep yourself awake — no moving about, talking, or singing. Just respond to the light.”',
      'instructions.pvt': '“A counter will appear at irregular intervals and start counting up. Press the space bar as fast as you can to stop it. Your reaction time is then shown. Do not press when there is no counter.”',

      'btn.start': 'Start task',
      'btn.startCal': 'Check display & start',
      'btn.begin': 'Begin trial',
      'btn.back': 'Back',
      'btn.again': 'New trial',
      'btn.end': 'End',
      'btn.exportRaw': 'Export raw trials (CSV)',
      'btn.exportPm': 'Export per-minute (CSV)',
      'btn.exportSummary': 'Export summary (CSV)',
      'btn.allRaw': 'All raw',
      'btn.allPm': 'All per-minute',
      'btn.allSummary': 'All summaries',

      'task.hintBsrt': 'Press the space bar, or tap, each time the light appears.',
      'task.hintPvt': 'Press the space bar as fast as you can to stop the counter.',
      'countdown.hint': 'Get comfortable. Respond to every stimulus.',

      'results.heading': 'Result',
      'results.sleepOnset': 'Sleep onset',
      'results.overview': 'Test overview',
      'results.totalTrials': 'Total trial run',
      'results.hitRatio': 'Hit ratio',
      'results.hits': 'Hits',
      'results.misses': 'Misses',
      'results.lapses': 'Lapses',
      'results.falseStarts': 'False starts',
      'results.ep12': 'Error profile EP1–2',
      'results.ep36': 'Error profile EP3–6',
      'results.ep7': 'Error profile EP7+',
      'results.rt': 'Total test — reaction time (ms)',
      'results.rs': 'Total test — reaction speed (1000/RT)',
      'results.raw': 'Raw',
      'results.corrected': 'Corrected',
      'results.average': 'Average',
      'results.median': 'Median',
      'results.sd': 'SD',
      'results.fastest': '10% fastest',
      'results.slowest': '10% slowest',
      'results.perMinute': 'Per minute',
      'results.integrity': 'Response integrity',
      'results.integrityOk': 'Response pattern looks normal.',
      'results.integrityFlag': 'Repeated or continuous pressing detected — review before using this trial.',
      'results.kss': 'Karolinska Sleepiness Scale',
      'results.kssBefore': 'Before',
      'results.kssAfter': 'After'
    },

    fr: {
      'app.subtitle': 'Tâche comportementale de résistance au sommeil et tâche de vigilance psychomotrice.',
      'lang.label': 'Langue',

      'participant.heading': 'Participant',
      'participant.id': 'Identifiant du participant',
      'participant.name': 'Nom',
      'participant.birth': 'Date de naissance',
      'participant.education': 'Niveau d’études',
      'participant.address': 'Adresse',
      'participant.session': 'Session / condition',
      'participant.trial': 'Numéro d’essai',
      'participant.privacy': 'Le nom, l’adresse et la date de naissance identifient directement la personne. Ces données sont conservées sur cet appareil et inscrites dans chaque fichier exporté. Pour la recherche, préférez un identifiant pseudonymisé et conservez la clé d’identification séparément, conformément à votre plan de protection des données. La date et l’heure du test sont enregistrées automatiquement.',

      'task.heading': 'Tâche',
      'task.mode': 'Mode',
      'task.mode.bsrt': 'BSRT / OSLER — intervalle fixe',
      'task.mode.pvt': 'PVT — intervalle variable',
      'task.isi': 'Intervalle entre stimuli (ms)',
      'task.isiSet': 'Intervalles (secondes, séparés par des virgules)',
      'task.block': 'Durée du bloc (s)',
      'task.seed': 'Graine du programme (vide = aléatoire)',

      'protocol.heading': 'Protocole',
      'protocol.stim': 'Durée du stimulus / fenêtre de réponse (ms)',
      'protocol.max': 'Durée maximale (min)',
      'protocol.lapse': 'Seuil de défaillance (ms)',
      'protocol.criterion': 'Omissions consécutives = endormissement',
      'protocol.criterionOn': 'Arrêter le test dès le critère d’endormissement',
      'protocol.alarm': 'Émettre une alarme lorsque le critère d’endormissement est atteint',
      'protocol.fullscreen': 'Plein écran pendant la tâche',
      'protocol.explain': 'Une réponse dans la fenêtre est une réussite ; une réussite lente (au-delà du seuil) est une défaillance. Une réponse plus tardive que la fenêtre est une omission — son temps de réaction brut est conservé, mais elle ne réinitialise pas la série d’omissions. Une réponse en deçà du seuil de faux départ est un faux départ.',

      'kss.heading': 'Échelle de somnolence de Karolinska',
      'kss.when': 'Administrer la KSS',
      'kss.when.none': 'Pas du tout',
      'kss.when.before': 'Avant la tâche',
      'kss.when.after': 'Après la tâche',
      'kss.when.both': 'Avant et après',
      'kss.question': 'À quel point vous sentez-vous somnolent en ce moment ?',
      'kss.instruction': 'Choisissez l’affirmation qui décrit le mieux votre état actuel.',
      'kss.beforeTitle': 'Avant la tâche',
      'kss.afterTitle': 'Après la tâche',
      'kss.continue': 'Continuer',

      'corrections.heading': 'Corrections',
      'corrections.type': 'Correction',
      'corrections.none': 'Aucune',
      'corrections.falseStarts': 'Retirer uniquement les faux départs',
      'corrections.outliers': 'Retirer uniquement les valeurs aberrantes',
      'corrections.both': 'Retirer les deux',
      'corrections.threshold': 'Seuil de faux départ (ms)',
      'corrections.sd': 'Seuil des valeurs aberrantes (ET)',

      'instructions.heading': 'Consignes à lire au participant',
      'instructions.bsrt': '« Une lumière rouge clignotera toutes les quelques secondes. Chaque fois que vous la voyez, appuyez sur la barre d’espace. Essayez de rester éveillé, mais ne faites rien d’autre pour vous tenir éveillé — ne bougez pas, ne parlez pas, ne chantez pas. Répondez simplement à la lumière. »',
      'instructions.pvt': '« Un compteur apparaîtra à intervalles irréguliers et commencera à défiler. Appuyez sur la barre d’espace le plus vite possible pour l’arrêter. Votre temps de réaction s’affichera ensuite. N’appuyez pas lorsqu’il n’y a pas de compteur. »',

      'btn.start': 'Démarrer la tâche',
      'btn.startCal': 'Vérifier l’écran et démarrer',
      'btn.begin': 'Commencer l’essai',
      'btn.back': 'Retour',
      'btn.again': 'Nouvel essai',
      'btn.end': 'Terminer',
      'btn.exportRaw': 'Exporter les essais bruts (CSV)',
      'btn.exportPm': 'Exporter par minute (CSV)',
      'btn.exportSummary': 'Exporter le résumé (CSV)',
      'btn.allRaw': 'Tout brut',
      'btn.allPm': 'Tout par minute',
      'btn.allSummary': 'Tous les résumés',

      'task.hintBsrt': 'Appuyez sur la barre d’espace, ou touchez l’écran, chaque fois que la lumière apparaît.',
      'task.hintPvt': 'Appuyez sur la barre d’espace le plus vite possible pour arrêter le compteur.',
      'countdown.hint': 'Installez-vous confortablement. Répondez à chaque stimulus.',

      'results.heading': 'Résultat',
      'results.sleepOnset': 'Endormissement',
      'results.overview': 'Aperçu du test',
      'results.totalTrials': 'Nombre total d’essais',
      'results.hitRatio': 'Taux de réussite',
      'results.hits': 'Réussites',
      'results.misses': 'Omissions',
      'results.lapses': 'Défaillances',
      'results.falseStarts': 'Faux départs',
      'results.ep12': 'Profil d’erreur EP1–2',
      'results.ep36': 'Profil d’erreur EP3–6',
      'results.ep7': 'Profil d’erreur EP7+',
      'results.rt': 'Test complet — temps de réaction (ms)',
      'results.rs': 'Test complet — vitesse de réaction (1000/TR)',
      'results.raw': 'Brut',
      'results.corrected': 'Corrigé',
      'results.average': 'Moyenne',
      'results.median': 'Médiane',
      'results.sd': 'ET',
      'results.fastest': '10 % les plus rapides',
      'results.slowest': '10 % les plus lents',
      'results.perMinute': 'Par minute',
      'results.integrity': 'Intégrité des réponses',
      'results.integrityOk': 'Le profil de réponse paraît normal.',
      'results.integrityFlag': 'Appuis répétés ou continus détectés — vérifiez avant d’utiliser cet essai.',
      'results.kss': 'Échelle de somnolence de Karolinska',
      'results.kssBefore': 'Avant',
      'results.kssAfter': 'Après'
    },

    nl: {
      'app.subtitle': 'Gedragsmatige slaapweerstandstaak en psychomotorische vigilantietaak.',
      'lang.label': 'Taal',

      'participant.heading': 'Deelnemer',
      'participant.id': 'Deelnemersnummer',
      'participant.name': 'Naam',
      'participant.birth': 'Geboortedatum',
      'participant.education': 'Opleidingsniveau',
      'participant.address': 'Adres',
      'participant.session': 'Sessie / conditie',
      'participant.trial': 'Testnummer',
      'participant.privacy': 'Naam, adres en geboortedatum zijn direct identificerend. Ze worden op dit apparaat bewaard en in elk geëxporteerd bestand geschreven. Gebruik voor onderzoek bij voorkeur een pseudoniem deelnemersnummer en bewaar de sleutel apart volgens uw gegevensbeschermingsplan. Datum en tijd van de test worden automatisch vastgelegd.',

      'task.heading': 'Taak',
      'task.mode': 'Modus',
      'task.mode.bsrt': 'BSRT / OSLER — vast interval',
      'task.mode.pvt': 'PVT — variabel interval',
      'task.isi': 'Interval tussen stimuli (ms)',
      'task.isiSet': 'Intervallen (seconden, komma-gescheiden)',
      'task.block': 'Bloklengte (s)',
      'task.seed': 'Seed van het schema (leeg = willekeurig)',

      'protocol.heading': 'Protocol',
      'protocol.stim': 'Stimulusduur / responsvenster (ms)',
      'protocol.max': 'Maximale duur (min)',
      'protocol.lapse': 'Drempel voor uitval (ms)',
      'protocol.criterion': 'Opeenvolgende missers = slaapaanvang',
      'protocol.criterionOn': 'Test vroegtijdig beëindigen bij het slaapaanvangscriterium',
      'protocol.alarm': 'Alarm laten klinken wanneer het slaapaanvangscriterium wordt bereikt',
      'protocol.fullscreen': 'Volledig scherm tijdens de taak',
      'protocol.explain': 'Een respons binnen het venster is een treffer; een trage treffer (voorbij de drempel) is een uitval. Een respons later dan het venster is een misser — de ruwe reactietijd wordt wel bewaard, maar de reeks missers wordt niet gereset. Een respons onder de drempel voor valse starts is een valse start.',

      'kss.heading': 'Karolinska-slaperigheidsschaal',
      'kss.when': 'KSS afnemen',
      'kss.when.none': 'Niet',
      'kss.when.before': 'Vóór de taak',
      'kss.when.after': 'Na de taak',
      'kss.when.both': 'Vóór en na',
      'kss.question': 'Hoe slaperig voelt u zich op dit moment?',
      'kss.instruction': 'Kies de uitspraak die het best beschrijft hoe u zich nu voelt.',
      'kss.beforeTitle': 'Vóór de taak',
      'kss.afterTitle': 'Na de taak',
      'kss.continue': 'Doorgaan',

      'corrections.heading': 'Correcties',
      'corrections.type': 'Correctie',
      'corrections.none': 'Geen',
      'corrections.falseStarts': 'Alleen valse starts verwijderen',
      'corrections.outliers': 'Alleen uitschieters verwijderen',
      'corrections.both': 'Beide verwijderen',
      'corrections.threshold': 'Drempel voor valse start (ms)',
      'corrections.sd': 'Afkapwaarde uitschieters (SD)',

      'instructions.heading': 'Instructies om aan de deelnemer voor te lezen',
      'instructions.bsrt': '“Om de paar seconden knippert een rood lampje. Druk telkens op de spatiebalk zodra u het ziet. Probeer wakker te blijven, maar doe verder niets om wakker te blijven — niet bewegen, praten of zingen. Reageer alleen op het lampje.”',
      'instructions.pvt': '“Op onregelmatige momenten verschijnt een teller die begint te lopen. Druk zo snel mogelijk op de spatiebalk om hem te stoppen. Daarna ziet u uw reactietijd. Druk niet wanneer er geen teller is.”',

      'btn.start': 'Taak starten',
      'btn.startCal': 'Scherm controleren en starten',
      'btn.begin': 'Test beginnen',
      'btn.back': 'Terug',
      'btn.again': 'Nieuwe test',
      'btn.end': 'Beëindigen',
      'btn.exportRaw': 'Ruwe gegevens exporteren (CSV)',
      'btn.exportPm': 'Per minuut exporteren (CSV)',
      'btn.exportSummary': 'Samenvatting exporteren (CSV)',
      'btn.allRaw': 'Alle ruwe',
      'btn.allPm': 'Alle per minuut',
      'btn.allSummary': 'Alle samenvattingen',

      'task.hintBsrt': 'Druk op de spatiebalk, of tik, telkens als het lampje verschijnt.',
      'task.hintPvt': 'Druk zo snel mogelijk op de spatiebalk om de teller te stoppen.',
      'countdown.hint': 'Ga comfortabel zitten. Reageer op elke stimulus.',

      'results.heading': 'Resultaat',
      'results.sleepOnset': 'Slaapaanvang',
      'results.overview': 'Overzicht van de test',
      'results.totalTrials': 'Totaal aantal trials',
      'results.hitRatio': 'Trefferpercentage',
      'results.hits': 'Treffers',
      'results.misses': 'Missers',
      'results.lapses': 'Uitvallen',
      'results.falseStarts': 'Valse starts',
      'results.ep12': 'Foutprofiel EP1–2',
      'results.ep36': 'Foutprofiel EP3–6',
      'results.ep7': 'Foutprofiel EP7+',
      'results.rt': 'Hele test — reactietijd (ms)',
      'results.rs': 'Hele test — reactiesnelheid (1000/RT)',
      'results.raw': 'Ruw',
      'results.corrected': 'Gecorrigeerd',
      'results.average': 'Gemiddelde',
      'results.median': 'Mediaan',
      'results.sd': 'SD',
      'results.fastest': '10% snelste',
      'results.slowest': '10% traagste',
      'results.perMinute': 'Per minuut',
      'results.integrity': 'Betrouwbaarheid van de responsen',
      'results.integrityOk': 'Het responspatroon lijkt normaal.',
      'results.integrityFlag': 'Herhaald of continu drukken vastgesteld — controleer voordat u deze test gebruikt.',
      'results.kss': 'Karolinska-slaperigheidsschaal',
      'results.kssBefore': 'Vóór',
      'results.kssAfter': 'Na'
    },

    de: {
      'app.subtitle': 'Verhaltensbasierter Schlafwiderstandstest und psychomotorischer Vigilanztest.',
      'lang.label': 'Sprache',

      'participant.heading': 'Teilnehmer',
      'participant.id': 'Teilnehmer-ID',
      'participant.name': 'Name',
      'participant.birth': 'Geburtsdatum',
      'participant.education': 'Bildungsabschluss',
      'participant.address': 'Adresse',
      'participant.session': 'Sitzung / Bedingung',
      'participant.trial': 'Durchgangsnummer',
      'participant.privacy': 'Name, Adresse und Geburtsdatum sind direkt identifizierend. Sie werden auf diesem Gerät gespeichert und in jede exportierte Datei geschrieben. Verwenden Sie für Forschungszwecke bevorzugt eine pseudonyme Teilnehmer-ID und bewahren Sie den Identitätsschlüssel gemäß Ihrem Datenschutzkonzept getrennt auf. Datum und Uhrzeit der Testung werden automatisch erfasst.',

      'task.heading': 'Aufgabe',
      'task.mode': 'Modus',
      'task.mode.bsrt': 'BSRT / OSLER — festes Intervall',
      'task.mode.pvt': 'PVT — variables Intervall',
      'task.isi': 'Interstimulusintervall (ms)',
      'task.isiSet': 'Intervalle (Sekunden, durch Komma getrennt)',
      'task.block': 'Blocklänge (s)',
      'task.seed': 'Startwert des Ablaufplans (leer = zufällig)',

      'protocol.heading': 'Protokoll',
      'protocol.stim': 'Stimulusdauer / Antwortfenster (ms)',
      'protocol.max': 'Maximale Dauer (Min.)',
      'protocol.lapse': 'Schwelle für Aussetzer (ms)',
      'protocol.criterion': 'Aufeinanderfolgende Auslassungen = Schlafbeginn',
      'protocol.criterionOn': 'Test beim Schlafbeginn-Kriterium vorzeitig beenden',
      'protocol.alarm': 'Alarm auslösen, wenn das Schlafbeginn-Kriterium erreicht wird',
      'protocol.fullscreen': 'Vollbild während der Aufgabe',
      'protocol.explain': 'Eine Antwort innerhalb des Fensters ist ein Treffer; ein langsamer Treffer (jenseits der Schwelle) ist ein Aussetzer. Eine Antwort nach dem Fenster ist eine Auslassung — ihre Rohreaktionszeit wird weiterhin erfasst, sie setzt die Serie von Auslassungen jedoch nicht zurück. Eine Antwort unterhalb der Fehlstart-Schwelle ist ein Fehlstart.',

      'kss.heading': 'Karolinska-Schläfrigkeitsskala',
      'kss.when': 'KSS erheben',
      'kss.when.none': 'Gar nicht',
      'kss.when.before': 'Vor der Aufgabe',
      'kss.when.after': 'Nach der Aufgabe',
      'kss.when.both': 'Vorher und nachher',
      'kss.question': 'Wie schläfrig fühlen Sie sich im Moment?',
      'kss.instruction': 'Wählen Sie die Aussage, die Ihren derzeitigen Zustand am besten beschreibt.',
      'kss.beforeTitle': 'Vor der Aufgabe',
      'kss.afterTitle': 'Nach der Aufgabe',
      'kss.continue': 'Weiter',

      'corrections.heading': 'Korrekturen',
      'corrections.type': 'Korrektur',
      'corrections.none': 'Keine',
      'corrections.falseStarts': 'Nur Fehlstarts entfernen',
      'corrections.outliers': 'Nur Ausreißer entfernen',
      'corrections.both': 'Beide entfernen',
      'corrections.threshold': 'Fehlstart-Schwelle (ms)',
      'corrections.sd': 'Ausreißergrenze (SD)',

      'instructions.heading': 'Anweisungen zum Vorlesen',
      'instructions.bsrt': '„Alle paar Sekunden blinkt ein rotes Licht auf. Drücken Sie jedes Mal die Leertaste, wenn Sie es sehen. Versuchen Sie wach zu bleiben, tun Sie aber sonst nichts, um sich wach zu halten — nicht bewegen, sprechen oder singen. Reagieren Sie nur auf das Licht.“',
      'instructions.pvt': '„In unregelmäßigen Abständen erscheint ein Zähler und beginnt zu laufen. Drücken Sie so schnell wie möglich die Leertaste, um ihn zu stoppen. Anschließend wird Ihre Reaktionszeit angezeigt. Drücken Sie nicht, wenn kein Zähler zu sehen ist.“',

      'btn.start': 'Aufgabe starten',
      'btn.startCal': 'Bildschirm prüfen und starten',
      'btn.begin': 'Durchgang beginnen',
      'btn.back': 'Zurück',
      'btn.again': 'Neuer Durchgang',
      'btn.end': 'Beenden',
      'btn.exportRaw': 'Rohdaten exportieren (CSV)',
      'btn.exportPm': 'Pro Minute exportieren (CSV)',
      'btn.exportSummary': 'Zusammenfassung exportieren (CSV)',
      'btn.allRaw': 'Alle Rohdaten',
      'btn.allPm': 'Alle pro Minute',
      'btn.allSummary': 'Alle Zusammenfassungen',

      'task.hintBsrt': 'Drücken Sie die Leertaste oder tippen Sie, sobald das Licht erscheint.',
      'task.hintPvt': 'Drücken Sie so schnell wie möglich die Leertaste, um den Zähler zu stoppen.',
      'countdown.hint': 'Machen Sie es sich bequem. Reagieren Sie auf jeden Reiz.',

      'results.heading': 'Ergebnis',
      'results.sleepOnset': 'Schlafbeginn',
      'results.overview': 'Testübersicht',
      'results.totalTrials': 'Durchgänge insgesamt',
      'results.hitRatio': 'Trefferquote',
      'results.hits': 'Treffer',
      'results.misses': 'Auslassungen',
      'results.lapses': 'Aussetzer',
      'results.falseStarts': 'Fehlstarts',
      'results.ep12': 'Fehlerprofil EP1–2',
      'results.ep36': 'Fehlerprofil EP3–6',
      'results.ep7': 'Fehlerprofil EP7+',
      'results.rt': 'Gesamttest — Reaktionszeit (ms)',
      'results.rs': 'Gesamttest — Reaktionsgeschwindigkeit (1000/RZ)',
      'results.raw': 'Roh',
      'results.corrected': 'Korrigiert',
      'results.average': 'Mittelwert',
      'results.median': 'Median',
      'results.sd': 'SD',
      'results.fastest': '10 % schnellste',
      'results.slowest': '10 % langsamste',
      'results.perMinute': 'Pro Minute',
      'results.integrity': 'Antwortintegrität',
      'results.integrityOk': 'Das Antwortmuster wirkt unauffällig.',
      'results.integrityFlag': 'Wiederholtes oder dauerhaftes Drücken erkannt — vor Verwendung dieses Durchgangs prüfen.',
      'results.kss': 'Karolinska-Schläfrigkeitsskala',
      'results.kssBefore': 'Vorher',
      'results.kssAfter': 'Nachher'
    }
  };

  var current = 'en';

  /* Falls back to English for any key a translation is missing, so an
   * incomplete translation degrades to English rather than to a blank. */
  function t(key) {
    var table = STRINGS[current] || STRINGS.en;
    if (table[key] !== undefined) return table[key];
    return STRINGS.en[key] !== undefined ? STRINGS.en[key] : key;
  }

  function kssAnchors(lang) {
    return KSS[lang || current] || KSS.en;
  }

  function setLanguage(code) {
    current = STRINGS[code] ? code : 'en';
    return current;
  }

  function getLanguage() { return current; }

  function applyTranslations(root) {
    var scope = root || document;
    scope.querySelectorAll('[data-i18n]').forEach(function (el) {
      el.textContent = t(el.getAttribute('data-i18n'));
    });
    scope.querySelectorAll('[data-i18n-ph]').forEach(function (el) {
      el.setAttribute('placeholder', t(el.getAttribute('data-i18n-ph')));
    });
    scope.querySelectorAll('[data-i18n-html]').forEach(function (el) {
      el.innerHTML = t(el.getAttribute('data-i18n-html'));
    });
  }

  return {
    LANGUAGES: LANGUAGES,
    t: t,
    kssAnchors: kssAnchors,
    setLanguage: setLanguage,
    getLanguage: getLanguage,
    applyTranslations: applyTranslations
  };
});
