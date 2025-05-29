// Initialize particles.js with the configuration directly
particlesJS('particles-js', {
    "particles": {
        "number": {
            "value": 80,
            "density": {
                "enable": true,
                "value_area": 800
            }
        },
        "color": {
            "value": "#00ccff"
        },
        "shape": {
            "type": "circle",
            "stroke": {
                "width": 0,
                "color": "#000000"
            }
        },
        "opacity": {
            "value": 0.5,
            "random": true,
            "anim": {
                "enable": false
            }
        },
        "size": {
            "value": 3,
            "random": true,
            "anim": {
                "enable": false
            }
        },
        "line_linked": {
            "enable": true,
            "distance": 150,
            "color": "#00ccff",
            "opacity": 0.4,
            "width": 1
        },
        "move": {
            "enable": true,
            "speed": 2,
            "direction": "none",
            "random": false,
            "straight": false,
            "out_mode": "out",
            "bounce": false
        }
    },
    "interactivity": {
        "detect_on": "canvas",
        "events": {
            "onhover": {
                "enable": true,
                "mode": "repulse"
            },
            "onclick": {
                "enable": true,
                "mode": "push"
            },
            "resize": true
        },
        "modes": {
            "repulse": {
                "distance": 100,
                "duration": 0.4
            },
            "push": {
                "particles_nb": 4
            }
        }
    },
    "retina_detect": true
}, function() {
    console.log('particles.js loaded');
});

// Typing Animation
document.addEventListener('DOMContentLoaded', function() {
    const typingText = document.getElementById('typing-text');
    const educationText = document.getElementById('education-text');
    const languagesText = document.getElementById('languages-text');

    const texts = [
        "Welcome to the YABS repository",
        "This repository is a hub of many engineering models, designed to be fast and practical and give common sense for engineering judgment.",
        "Use them to get some ideas about how things work. Explore the matrix of features and feel free to contribute.",
        "And remember it is a work in progress!",
    ];

    let currentTextIndex = 0;
    let currentCharIndex = 0;
    let isTyping = false;

    function type() {
        if (currentTextIndex >= texts.length) return; // Stop if all texts are typed

        isTyping = true;
        const currentText = texts[currentTextIndex];
        let targetElement;

        // Determine which element to type into based on the current text index
        if (currentTextIndex === 0 || currentTextIndex === 1) {
            targetElement = typingText;
        } else if (currentTextIndex === 2) {
            targetElement = educationText;
        } else if (currentTextIndex === 3) {
            targetElement = languagesText;
        }

        // Add the next character
        if (currentCharIndex < currentText.length) {
            targetElement.textContent = currentText.substring(0, currentCharIndex + 1);
            currentCharIndex++;
            setTimeout(type, 50); // Adjust typing speed (50ms per character)
        } else {
            // Move to the next text after a short pause
            currentTextIndex++;
            currentCharIndex = 0;
            isTyping = false;
            setTimeout(type, 1000); // Pause for 1 second before typing the next sentence
        }
    }

    // Start typing when the page loads
    type();
});

// Toggle game visibility
function toggleGame() {
    const gameContainer = document.getElementById('gameContainer');
    const gameCanvas = document.getElementById('gameCanvas');
    if (gameContainer.style.display === 'none') {
        gameContainer.style.display = 'block';
        startGame(gameCanvas); // Start the game when opened
    } else {
        gameContainer.style.display = 'none';
        // Reset game state if needed
        resetGame();
    }
}

// Toggle contact form visibility
function toggleContactForm() {
    const contactFormContainer = document.getElementById('contactFormContainer');
    if (contactFormContainer.style.display === 'none') {
        contactFormContainer.style.display = 'flex';
    } else {
        contactFormContainer.style.display = 'none';
        // Reset the form when closing
        document.getElementById('contactForm').reset();
    }
}
